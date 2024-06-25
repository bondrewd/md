const std = @import("std");
const F = @import("precision.zig").F;
const Cfg = @import("cfg.zig").Cfg;
const Top = @import("top.zig").Top;
const Par = @import("par.zig").Par;
const Particle = @import("particle.zig").Particle;
const Class = @import("par.zig").Class;
const LJ = @import("par.zig").LJ;
const Crd = @import("crd.zig").Crd;
const BOLTZMANN = @import("constants.zig").BOLTZMANN;

const Energy = struct {
    kinetic: F,
    lj: F,
};

const Side = struct {
    full: F,
    half: F,
};

const Boundary = struct {
    x: ?Side = null,
    y: ?Side = null,
    z: ?Side = null,

    const Self = @This();

    fn wrap(self: *Self, r: *[3]F) void {
        if (self.x) |x| {
            if (r[0] > x.half or r[0] < -x.half) r[0] -= x.full * @round(r[0] / x.full);
        }
        if (self.y) |y| {
            if (r[1] > y.half or r[1] < -y.half) r[1] -= y.full * @round(r[1] / y.full);
        }
        if (self.z) |z| {
            if (r[2] > z.half or r[2] < -z.half) r[2] -= z.full * @round(r[2] / z.full);
        }
    }
};

pub const System = struct {
    n: u64,
    r: [][3]F,
    v: [][3]F,
    f: [][3]F,
    m: []F,
    q: []F,
    e: []F,
    s: []F,
    id: []u64,
    energy: Energy,
    boundary: Boundary,
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, n: u64) !Self {
        const new_system = Self{
            .n = n,
            .r = try allocator.alloc([3]F, n),
            .v = try allocator.alloc([3]F, n),
            .f = try allocator.alloc([3]F, n),
            .m = try allocator.alloc(F, n),
            .q = try allocator.alloc(F, n),
            .e = try allocator.alloc(F, n),
            .s = try allocator.alloc(F, n),
            .id = try allocator.alloc(u64, n),
            .energy = .{ .kinetic = undefined, .lj = undefined },
            .boundary = .{},
            .allocator = allocator,
        };
        return new_system;
    }

    pub fn deinit(self: *System) void {
        self.allocator.free(self.r);
        self.allocator.free(self.v);
        self.allocator.free(self.f);
        self.allocator.free(self.m);
        self.allocator.free(self.q);
        self.allocator.free(self.e);
        self.allocator.free(self.s);
        self.allocator.free(self.id);
    }

    pub fn setFromCfg(self: *Self, cfg: Cfg) void {
        if (cfg.boundary) |boundary| {
            if (boundary.x) |x| self.boundary.x = .{ .full = x, .half = x / 2 };
            if (boundary.y) |y| self.boundary.y = .{ .full = y, .half = y / 2 };
            if (boundary.z) |z| self.boundary.z = .{ .full = z, .half = z / 2 };
        }
    }

    pub fn setFromCrd(self: *Self, crd: Crd) void {
        for (crd.coordinates, self.r, self.id) |coordinate, *r, *id| {
            r.* = .{ coordinate.x, coordinate.y, coordinate.z };
            id.* = coordinate.id;
        }
    }

    pub fn setFromTopPar(self: *Self, top: Top, par: Par) void {
        for (self.id, self.m, self.q, self.e, self.s) |id, *m, *q, *e, *s| {
            var particle: ?Particle = null;
            for (top.particles) |top_particle| {
                if (top_particle.id == id) {
                    particle = top_particle;
                    break;
                }
            }

            if (particle == null) {
                std.debug.print("[ERROR] Particle with id {d} not found\n", .{id});
                std.process.exit(0);
            }

            var class: ?Class = null;
            for (par.classes) |par_class| {
                if (std.mem.eql(u8, par_class.name, particle.?.name)) {
                    class = par_class;
                    break;
                }
            }

            if (class == null) {
                std.debug.print("[ERROR] Class with name {s} not found\n", .{particle.?.name});
                std.process.exit(0);
            }

            m.* = class.?.mass;
            q.* = class.?.charge;

            var lj: ?LJ = null;
            for (par.lj) |par_lj| {
                if (std.mem.eql(u8, par_lj.name, particle.?.name)) {
                    lj = par_lj;
                    break;
                }
            }

            if (lj == null) {
                std.debug.print("[ERROR] lj parameters for name {s} not found\n", .{particle.?.name});
                std.process.exit(0);
            }

            e.* = lj.?.epsilon;
            s.* = lj.?.sigma;
        }
    }

    pub fn setRandomVelocities(self: *Self, rng: std.Random) void {
        for (self.v) |*v| {
            v.* = .{
                2 * rng.float(F) - 1,
                2 * rng.float(F) - 1,
                2 * rng.float(F) - 1,
            };
        }
    }

    pub fn setRandomVelocitiesWithTemperature(self: *Self, rng: std.Random, temperature: F) void {
        self.setRandomVelocities(rng);
        const target_temperature = temperature;
        const current_temperature = self.measureTemperature();
        const scaling_factor = std.math.sqrt(target_temperature / current_temperature);
        for (self.v) |*v| {
            v[0] *= scaling_factor;
            v[1] *= scaling_factor;
            v[2] *= scaling_factor;
        }
    }

    pub fn zeroEnergy(self: *Self) void {
        self.energy.kinetic = 0.0;
        self.energy.lj = 0.0;
    }

    pub fn zeroForce(self: *Self) void {
        for (self.f) |*f| f.* = .{ 0, 0, 0 };
    }

    pub fn wrapCoordinates(self: *Self) void {
        for (self.r) |*r| self.boundary.wrap(r);
    }

    pub fn updateForceLJ(self: *Self) void {
        for (0..self.n - 1) |i| {
            const ri = self.r[i];
            const ei = self.e[i];
            const si = self.s[i];
            for (i + 1..self.n) |j| {
                const ej = self.e[j];
                const sj = self.s[j];
                const e = @sqrt(ei * ej);
                const s = (si + sj) / 2;

                const rj = self.r[j];
                var dr = [3]F{ rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2] };
                self.boundary.wrap(&dr);

                const s2 = s * s;
                const c2 = s2 / (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
                const c4 = c2 * c2;
                const c6 = c4 * c2;
                const c8 = c4 * c4;
                const force = 48 * e * c8 * (c6 - 0.5) / s2;

                self.f[i][0] -= force * dr[0];
                self.f[i][1] -= force * dr[1];
                self.f[i][2] -= force * dr[2];

                self.f[j][0] += force * dr[0];
                self.f[j][1] += force * dr[1];
                self.f[j][2] += force * dr[2];
            }
        }
    }

    pub fn updateEnergyLJ(self: *Self) void {
        for (0..self.n - 1) |i| {
            const ri = self.r[i];
            const ei = self.e[i];
            const si = self.s[i];
            for (i + 1..self.n) |j| {
                const ej = self.e[j];
                const sj = self.s[j];
                const e = @sqrt(ei * ej);
                const s = (si + sj) / 2;

                const rj = self.r[j];
                var dr = [3]F{ rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2] };
                self.boundary.wrap(&dr);

                const s2 = s * s;
                const c2 = s2 / (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
                const c4 = c2 * c2;
                const c6 = c4 * c2;
                const energy = 4 * e * c6 * (c6 - 1);

                self.energy.lj += energy;
            }
        }
    }

    pub fn updateEnergyKinetic(self: *Self) void {
        for (self.v, self.m) |v, m| {
            self.energy.kinetic += (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * m / 2;
        }
        self.energy.kinetic /= 418.4;
    }

    pub fn measureTemperature(self: *const Self) F {
        var sum: F = 0.0;
        for (self.v, self.m) |v, m| sum += (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * m / (3.0 * BOLTZMANN);
        sum /= 418.4;
        const temperature = sum / @as(F, @floatFromInt(self.n));
        return temperature;
    }
};
