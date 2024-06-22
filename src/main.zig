const std = @import("std");

const BOLTZMANN = 1.987204259e-3;

const InputBlock = struct {
    crd: []u8,
    par: []u8,
    top: []u8,

    const Self = @This();

    fn log(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] INPUT\n", .{});
        try writer.print("[INFO] crd = {s}\n", .{self.crd});
        try writer.print("[INFO] par = {s}\n", .{self.par});
        try writer.print("[INFO] top = {s}\n", .{self.top});
    }
};

const OutputBlock = struct {
    crd: ?[]u8 = null,
    xyz: ?[]u8 = null,

    const Self = @This();

    fn log(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] OUTPUT\n", .{});
        if (self.crd) |crd| {
            try writer.print("[INFO] crd = {s}\n", .{crd});
        } else {
            try writer.print("[INFO] crd =\n", .{});
        }
        if (self.xyz) |xyz| {
            try writer.print("[INFO] xyz = {s}\n", .{xyz});
        } else {
            try writer.print("[INFO] xyz =\n", .{});
        }
    }
};

const EnsembleKind = enum {
    nve,
    nvt,
    npt,
};

const ThermostatKind = enum {
    berendsen,
    langevin,
};

const BarostatKind = enum {
    berendsen,
    langevin,
};

const DynamicsBlock = struct {
    dt: f64,
    steps: u64,
    ensemble: ?EnsembleKind = null,
    temperature: ?f64 = null,
    pressure: ?f64 = null,
    thermostat: ?ThermostatKind = null,
    barostat: ?BarostatKind = null,

    const Self = @This();

    fn log(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] DYNAMICS\n", .{});
        try writer.print("[INFO] dt = {d}\n", .{self.dt});
        try writer.print("[INFO] steps = {d}\n", .{self.steps});
        if (self.ensemble) |ensemble| {
            try writer.print("[INFO] ensemble = {s}\n", .{@tagName(ensemble)});
        } else {
            try writer.print("[INFO] ensemble =\n", .{});
        }
        if (self.temperature) |temperature| {
            try writer.print("[INFO] temperature = {d}\n", .{temperature});
        } else {
            try writer.print("[INFO] temperature =\n", .{});
        }
        if (self.pressure) |pressure| {
            try writer.print("[INFO] pressure = {d}\n", .{pressure});
        } else {
            try writer.print("[INFO] pressure =\n", .{});
        }
        if (self.thermostat) |thermostat| {
            try writer.print("[INFO] thermostat = {s}\n", .{@tagName(thermostat)});
        } else {
            try writer.print("[INFO] thermostat =\n", .{});
        }
        if (self.barostat) |barostat| {
            try writer.print("[INFO] barostat = {s}\n", .{@tagName(barostat)});
        } else {
            try writer.print("[INFO] barostat =\n", .{});
        }
    }
};

const Rng = struct {
    seed: u64,

    const Self = @This();

    fn log(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] RNG\n", .{});
        try writer.print("[INFO] seed = {d}\n", .{self.seed});
    }
};

const BoundaryBlock = struct {
    x: ?f64 = null,
    y: ?f64 = null,
    z: ?f64 = null,

    const Self = @This();

    fn log(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] BOUNDARY\n", .{});
        if (self.x) |x| {
            try writer.print("[INFO] x = {d}\n", .{x});
        }
        if (self.y) |y| {
            try writer.print("[INFO] y = {d}\n", .{y});
        }
        if (self.z) |z| {
            try writer.print("[INFO] z = {d}\n", .{z});
        }
    }
};

const Cfg = struct {
    input: InputBlock,
    output: ?OutputBlock = null,
    dynamics: ?DynamicsBlock = null,
    boundary: ?BoundaryBlock = null,
    rng: ?Rng = null,

    const Self = @This();

    fn log(self: *const Self, writer: std.fs.File.Writer) !void {
        try self.input.log(writer);
        try writer.print("[INFO] \n", .{});
        if (self.output) |output| {
            try output.log(writer);
            try writer.print("[INFO] \n", .{});
        }
        if (self.dynamics) |dynamics| {
            try dynamics.log(writer);
            try writer.print("[INFO] \n", .{});
        }
        if (self.rng) |rng| {
            try rng.log(writer);
            try writer.print("[INFO] \n", .{});
        }
        if (self.boundary) |boundary| {
            try boundary.log(writer);
            try writer.print("[INFO] \n", .{});
        }
    }
};

const Coordinate = struct {
    id: u64,
    x: f64,
    y: f64,
    z: f64,
};

const Crd = struct {
    coordinates: []Coordinate,

    const Self = @This();

    fn log(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] CRD\n", .{});
        for (self.coordinates) |coordinate| {
            try writer.print("[INFO] id = {d}, x = {d}, y = {d}, z = {d}\n", .{ coordinate.id, coordinate.x, coordinate.y, coordinate.z });
        }
        try writer.print("[INFO]\n", .{});
    }
};

const Class = struct {
    name: []const u8,
    mass: f64,
    charge: f64,
};

const LJ = struct {
    name: []const u8,
    epsilon: f64,
    sigma: f64,
};

const Par = struct {
    classes: []const Class,
    lj: []const LJ,

    const Self = @This();

    fn log(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] CLASSES\n", .{});
        for (self.classes) |class| {
            try writer.print("[INFO] name = {s}, mass = {d}, charge = {d}\n", .{ class.name, class.mass, class.charge });
        }
        try writer.print("[INFO]\n", .{});

        try writer.print("[INFO] LJ\n", .{});
        for (self.lj) |lj| {
            try writer.print("[INFO] name = {s}, epsilon = {d}, sigma = {d}\n", .{ lj.name, lj.epsilon, lj.sigma });
        }
        try writer.print("[INFO]\n", .{});
    }
};

const Particle = struct {
    id: u64,
    name: []const u8,
};

const Top = struct {
    particles: []Particle,

    const Self = @This();

    fn log(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] TOP\n", .{});
        for (self.particles) |particle| {
            try writer.print("[INFO] id = {d}, class = {s}\n", .{ particle.id, particle.name });
        }
        try writer.print("[INFO]\n", .{});
    }
};

const Energy = struct {
    kinetic: f64,
    lj: f64,
};

const Side = struct {
    full: f64,
    half: f64,
};

const Boundary = struct {
    x: ?Side = null,
    y: ?Side = null,
    z: ?Side = null,

    const Self = @This();

    fn wrap(self: *Self, r: *[3]f64) void {
        if (self.x) |x| {
            if (r[0] > x.half or r[0] < -x.half) r[0] -= x.full * @round(r[0] / x.full);
        }
        if (self.y) |y| {
            if (r[0] > y.half or r[0] < -y.half) r[0] -= y.full * @round(r[0] / y.full);
        }
        if (self.z) |z| {
            if (r[0] > z.half or r[0] < -z.half) r[0] -= z.full * @round(r[0] / z.full);
        }
    }
};

const System = struct {
    n: u64,
    r: [][3]f64,
    v: [][3]f64,
    f: [][3]f64,
    m: []f64,
    q: []f64,
    e: []f64,
    s: []f64,
    id: []u64,
    energy: Energy,
    boundary: Boundary,
    allocator: std.mem.Allocator,

    const Self = @This();

    fn init(allocator: std.mem.Allocator, n: u64) !Self {
        const new_system = Self{
            .n = n,
            .r = try allocator.alloc([3]f64, n),
            .v = try allocator.alloc([3]f64, n),
            .f = try allocator.alloc([3]f64, n),
            .m = try allocator.alloc(f64, n),
            .q = try allocator.alloc(f64, n),
            .e = try allocator.alloc(f64, n),
            .s = try allocator.alloc(f64, n),
            .id = try allocator.alloc(u64, n),
            .energy = .{ .kinetic = undefined, .lj = undefined },
            .boundary = .{},
            .allocator = allocator,
        };
        return new_system;
    }

    fn deinit(self: *System) void {
        self.allocator.free(self.r);
        self.allocator.free(self.v);
        self.allocator.free(self.f);
        self.allocator.free(self.m);
        self.allocator.free(self.q);
        self.allocator.free(self.e);
        self.allocator.free(self.s);
        self.allocator.free(self.id);
    }

    fn setFromCfg(self: *Self, cfg: Cfg) void {
        if (cfg.boundary) |boundary| {
            if (boundary.x) |x| self.boundary.x = .{ .full = x, .half = x / 2 };
            if (boundary.y) |y| self.boundary.y = .{ .full = y, .half = y / 2 };
            if (boundary.z) |z| self.boundary.z = .{ .full = z, .half = z / 2 };
        }
    }

    fn setFromCrd(self: *Self, crd: Crd) void {
        for (crd.coordinates, self.r, self.id) |coordinate, *r, *id| {
            r.* = .{ coordinate.x, coordinate.y, coordinate.z };
            id.* = coordinate.id;
        }
    }

    fn setFromTopPar(self: *Self, top: Top, par: Par) void {
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

    fn setRandomVelocities(self: *Self, rng: std.Random) void {
        for (self.v) |*v| {
            v.* = .{
                2 * rng.float(f64) - 1,
                2 * rng.float(f64) - 1,
                2 * rng.float(f64) - 1,
            };
        }
    }

    fn setRandomVelocitiesWithTemperature(self: *Self, rng: std.Random, temperature: f64) void {
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

    fn zeroEnergy(self: *Self) void {
        self.energy.kinetic = 0.0;
        self.energy.lj = 0.0;
    }

    fn zeroForce(self: *Self) void {
        for (self.f) |*f| f.* = .{ 0, 0, 0 };
    }

    fn wrapCoordinates(self: *Self) void {
        for (self.r) |*r| self.boundary.wrap(r);
    }

    fn updateForceLJ(self: *Self) void {
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
                var dr = [3]f64{ rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2] };
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

    fn updateEnergyLJ(self: *Self) void {
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
                var dr = [3]f64{ rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2] };
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

    fn updateEnergyKinetic(self: *Self) void {
        for (self.v, self.m) |v, m| {
            self.energy.kinetic += (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * m / 2;
        }
        self.energy.kinetic /= 418.4;
    }

    fn measureTemperature(self: *const Self) f64 {
        var sum: f64 = 0.0;
        for (self.v, self.m) |v, m| sum += (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * m / (3.0 * BOLTZMANN);
        sum /= 418.4;
        const temperature = sum / @as(f64, @floatFromInt(self.n));
        return temperature;
    }

    fn log_dynamics_header(_: *Self, writer: std.fs.File.Writer) !void {
        try writer.print("[SIMU] {s:>10} {s:>10} {s:>10} {s:>10} {s:>10} {s:>10}\n", .{
            "STEP",
            "TEMP",
            "TOTAL",
            "POTENTIAL",
            "KINETIC",
            "LJ",
        });
    }

    fn log_dynamics(self: *Self, writer: std.fs.File.Writer, step: usize) !void {
        try writer.print("[SIMU] {d:>10} {d:>10.2} {d:>10.4} {d:>10.4} {d:>10.4} {d:>10.4}\n", .{
            step,
            self.measureTemperature(),
            self.energy.lj + self.energy.kinetic,
            self.energy.lj,
            self.energy.kinetic,
            self.energy.lj,
        });
    }
};

const XYZFile = struct {
    file: ?std.fs.File = null,
    writer: ?std.fs.File.Writer = null,

    const Self = @This();

    fn init(path: ?[]const u8) Self {
        var xyz_file = Self{};
        if (path) |xyz_file_path| {
            xyz_file.file = std.fs.cwd().createFile(xyz_file_path, .{ .truncate = true }) catch {
                std.debug.print("[ERROR] Unable to open file '{s}'\n", .{xyz_file_path});
                std.process.exit(0);
            };
            xyz_file.writer = xyz_file.file.?.writer();
        }
        return xyz_file;
    }

    fn deinit(self: *const Self) void {
        if (self.file) |file| {
            file.close();
        }
    }

    fn write_frame(self: *const Self, system: *const System) !void {
        if (self.writer) |writer| {
            try writer.print("{d}\n\n", .{system.n});
            for (system.r) |r| {
                try writer.print("Ar {d:8.3} {d:8.3} {d:8.3}\n", .{ r[0], r[1], r[2] });
            }
        }
    }
};

pub fn main() !void {
    const stdout = std.io.getStdOut().writer();
    const stderr = std.io.getStdErr().writer();

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        const deinit_status = gpa.deinit();
        if (deinit_status == .leak) std.testing.expect(false) catch @panic("TEST FAIL");
    }
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len != 2) {
        std.debug.print("[ERROR] Invalid argument\n", .{});
        std.debug.print("Usage: ./md CFG\n", .{});
        return;
    }

    // Read configuration file
    const cfg_file_path = args[1];
    const cfg_file = std.fs.cwd().openFile(cfg_file_path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("[ERROR] File '{s}' not found\n", .{cfg_file_path});
            return;
        },
        else => {
            std.debug.print("[ERROR] Unable to open file '{s}'\n", .{cfg_file_path});
            return;
        },
    };
    defer cfg_file.close();

    const cfg_data = try cfg_file.readToEndAlloc(allocator, 1 << 10);
    defer allocator.free(cfg_data);

    const cfg_parsed = try std.json.parseFromSlice(Cfg, allocator, cfg_data, .{});
    defer cfg_parsed.deinit();

    const cfg = cfg_parsed.value;
    try cfg.log(stderr);

    // Read coordinate file
    var crd_file_path: []u8 = undefined;
    if (std.fs.path.dirname(cfg_file_path)) |dir| {
        crd_file_path = try std.fs.path.join(allocator, &[_][]const u8{ dir, cfg.input.crd });
    } else {
        crd_file_path = try std.fs.path.join(allocator, &[_][]const u8{ ".", cfg.input.crd });
    }
    defer allocator.free(crd_file_path);
    const crd_file = std.fs.cwd().openFile(crd_file_path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("[ERROR] File '{s}' not found\n", .{crd_file_path});
            return;
        },
        else => {
            std.debug.print("[ERROR] Unable to open file '{s}'\n", .{crd_file_path});
            return;
        },
    };
    defer crd_file.close();

    const crd_data = try crd_file.readToEndAlloc(allocator, 1 << 20);
    defer allocator.free(crd_data);

    const crd_parsed = try std.json.parseFromSlice(Crd, allocator, crd_data, .{});
    defer crd_parsed.deinit();

    const crd = crd_parsed.value;
    try crd.log(stderr);

    // Read parameter file
    var par_file_path: []u8 = undefined;
    if (std.fs.path.dirname(cfg_file_path)) |dir| {
        par_file_path = try std.fs.path.join(allocator, &[_][]const u8{ dir, cfg.input.par });
    } else {
        par_file_path = try std.fs.path.join(allocator, &[_][]const u8{ ".", cfg.input.par });
    }
    defer allocator.free(par_file_path);
    const par_file = std.fs.cwd().openFile(par_file_path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("[ERROR] File '{s}' not found\n", .{par_file_path});
            return;
        },
        else => {
            std.debug.print("[ERROR] Unable to open file '{s}'\n", .{par_file_path});
            return;
        },
    };
    defer par_file.close();

    const par_data = try par_file.readToEndAlloc(allocator, 1 << 20);
    defer allocator.free(par_data);

    const par_parsed = try std.json.parseFromSlice(Par, allocator, par_data, .{});
    defer par_parsed.deinit();

    const par = par_parsed.value;
    try par.log(stderr);

    // Read parameter file
    var top_file_path: []u8 = undefined;
    if (std.fs.path.dirname(cfg_file_path)) |dir| {
        top_file_path = try std.fs.path.join(allocator, &[_][]const u8{ dir, cfg.input.top });
    } else {
        top_file_path = try std.fs.path.join(allocator, &[_][]const u8{ ".", cfg.input.top });
    }
    defer allocator.free(top_file_path);
    const top_file = std.fs.cwd().openFile(top_file_path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("[ERROR] File '{s}' not found\n", .{top_file_path});
            return;
        },
        else => {
            std.debug.print("[ERROR] Unable to open file '{s}'\n", .{top_file_path});
            return;
        },
    };
    defer top_file.close();

    const top_data = try top_file.readToEndAlloc(allocator, 1 << 20);
    defer allocator.free(top_data);

    const top_topsed = try std.json.parseFromSlice(Top, allocator, top_data, .{});
    defer top_topsed.deinit();

    const top = top_topsed.value;
    try top.log(stderr);

    // Initialize rng seed
    const seed = blk: {
        if (cfg.rng) |rng| {
            break :blk rng.seed;
        } else {
            var seed: u64 = undefined;
            try std.posix.getrandom(std.mem.asBytes(&seed));
            break :blk seed;
        }
    };
    var prng = std.rand.DefaultPrng.init(seed);
    const rng = prng.random();
    std.debug.print("[INFO] SETUP RNG\n", .{});
    std.debug.print("[INFO] seed = {d}\n", .{seed});
    std.debug.print("[INFO] type = {s}\n", .{@typeName(@TypeOf(prng))});
    std.debug.print("[INFO]\n", .{});

    // Initialize system
    var system = try System.init(allocator, crd.coordinates.len);
    defer system.deinit();
    system.setFromCfg(cfg);
    system.setFromCrd(crd);
    system.setFromTopPar(top, par);
    system.wrapCoordinates();

    // Initialize velocities
    if (cfg.dynamics) |dynamics| {
        if (dynamics.temperature) |temperature| {
            system.setRandomVelocitiesWithTemperature(rng, temperature);
        } else {
            system.setRandomVelocities(rng);
        }
    }

    // Initialize force
    system.zeroForce();
    system.updateForceLJ();

    // Initialize energy
    system.zeroEnergy();
    system.updateEnergyLJ();
    system.updateEnergyKinetic();

    std.debug.print("[INFO] SETUP SYSTEM\n", .{});
    std.debug.print("[INFO] number of particles: {d}\n", .{system.n});
    std.debug.print("[INFO] initial temperature: {d:.2} K\n", .{system.measureTemperature()});
    std.debug.print("[INFO] initial energy lj: {d:.2} kcal/mol\n", .{system.energy.lj});
    std.debug.print("[INFO] initial energy kinetic: {d:.2} kcal/mol\n", .{system.energy.kinetic});
    std.debug.print("[INFO] initial energy total: {d:.2} kcal/mol\n", .{system.energy.lj + system.energy.kinetic});
    std.debug.print("[INFO]\n", .{});

    // Run dynamics
    if (cfg.dynamics) |dynamics| {
        // Open output files
        var xyz_file_path: ?[]u8 = null;
        if (cfg.output) |output| {
            if (output.xyz) |xyz| {
                if (std.fs.path.dirname(cfg_file_path)) |dir| {
                    xyz_file_path = try std.fs.path.join(allocator, &[_][]const u8{ dir, xyz });
                } else {
                    xyz_file_path = try std.fs.path.join(allocator, &[_][]const u8{ ".", xyz });
                }
            }
        }
        defer if (xyz_file_path) |path| allocator.free(path);
        const xyz_file = XYZFile.init(xyz_file_path);
        defer xyz_file.deinit();
        try xyz_file.write_frame(&system);
        // Log header and first step
        try system.log_dynamics_header(stdout);
        try system.log_dynamics(stdout, 0);
        // Initialize variables
        const steps = dynamics.steps;
        const dt = dynamics.dt;
        for (1..steps + 1) |step| {
            // Update coordinates and velocities
            for (system.r, system.v, system.f, system.m) |*r, *v, f, m| {
                r[0] += v[0] * dt + 418.4 * f[0] * dt * dt / (2 * m);
                r[1] += v[1] * dt + 418.4 * f[1] * dt * dt / (2 * m);
                r[2] += v[2] * dt + 418.4 * f[2] * dt * dt / (2 * m);
                v[0] += 418.4 * f[0] * dt / (2 * m);
                v[1] += 418.4 * f[1] * dt / (2 * m);
                v[2] += 418.4 * f[2] * dt / (2 * m);
            }
            // Wrap coordinates
            system.wrapCoordinates();
            // Update forces
            system.zeroForce();
            system.updateForceLJ();
            // Update velocities
            for (system.v, system.f, system.m) |*v, f, m| {
                v[0] += 418.4 * f[0] * dt / (2 * m);
                v[1] += 418.4 * f[1] * dt / (2 * m);
                v[2] += 418.4 * f[2] * dt / (2 * m);
            }
            // Update energy
            system.zeroEnergy();
            system.updateEnergyLJ();
            system.updateEnergyKinetic();

            if (step % 100 == 0) try system.log_dynamics(stdout, step);
            if (step % 100 == 0) try xyz_file.write_frame(&system);
        }
    }
}
