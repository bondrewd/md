const std = @import("std");

const InputBlock = struct {
    crd: []u8,
    par: []u8,
    top: []u8,
};

const OutputBlock = struct {
    crd: ?[]u8 = null,
    trj: ?[]u8 = null,
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
    dt: ?f64 = null,
    steps: ?u64 = null,
    ensemble: ?EnsembleKind = null,
    temperature: ?f64 = null,
    pressure: ?f64 = null,
    thermostat: ?ThermostatKind = null,
    barostat: ?BarostatKind = null,
};

const BoundaryBlock = struct {
    x: ?f64 = null,
    y: ?f64 = null,
    z: ?f64 = null,
};

const Cfg = struct {
    input: InputBlock,
    output: ?OutputBlock = null,
    dynamics: ?DynamicsBlock = null,
    boundary: ?BoundaryBlock = null,

    const Self = @This();

    fn display(self: *const Self) void {
        std.debug.print("[INFO] INPUT\n", .{});
        std.debug.print("[INFO] crd = {s}\n", .{self.input.crd});
        std.debug.print("[INFO] par = {s}\n", .{self.input.par});
        std.debug.print("[INFO] top = {s}\n", .{self.input.top});
        std.debug.print("[INFO]\n", .{});

        if (self.output) |output| {
            std.debug.print("[INFO] OUTPUT\n", .{});
            if (output.crd) |crd| {
                std.debug.print("[INFO] crd = {s}\n", .{crd});
            } else {
                std.debug.print("[INFO] crd =\n", .{});
            }
            if (output.trj) |trj| {
                std.debug.print("[INFO] trj = {s}\n", .{trj});
            } else {
                std.debug.print("[INFO] trj =\n", .{});
            }
            std.debug.print("[INFO]\n", .{});
        }

        if (self.dynamics) |dynamics| {
            std.debug.print("[INFO] DYNAMICS\n", .{});
            if (dynamics.dt) |dt| {
                std.debug.print("[INFO] dt = {d}\n", .{dt});
            } else {
                std.debug.print("[INFO] dt =\n", .{});
            }
            if (dynamics.steps) |steps| {
                std.debug.print("[INFO] steps = {d}\n", .{steps});
            } else {
                std.debug.print("[INFO] steps =\n", .{});
            }
            if (dynamics.ensemble) |ensemble| {
                std.debug.print("[INFO] ensemble = {s}\n", .{@tagName(ensemble)});
            } else {
                std.debug.print("[INFO] ensemble =\n", .{});
            }
            if (dynamics.temperature) |temperature| {
                std.debug.print("[INFO] temperature = {d}\n", .{temperature});
            } else {
                std.debug.print("[INFO] temperature =\n", .{});
            }
            if (dynamics.pressure) |pressure| {
                std.debug.print("[INFO] pressure = {d}\n", .{pressure});
            } else {
                std.debug.print("[INFO] pressure =\n", .{});
            }
            if (dynamics.thermostat) |thermostat| {
                std.debug.print("[INFO] thermostat = {s}\n", .{@tagName(thermostat)});
            } else {
                std.debug.print("[INFO] thermostat =\n", .{});
            }
            if (dynamics.barostat) |barostat| {
                std.debug.print("[INFO] barostat = {s}\n", .{@tagName(barostat)});
            } else {
                std.debug.print("[INFO] barostat =\n", .{});
            }
            std.debug.print("[INFO]\n", .{});
        }

        if (self.boundary) |boundary| {
            std.debug.print("[INFO] BOUNDARY\n", .{});
            if (boundary.x) |x| {
                std.debug.print("[INFO] x = {d}\n", .{x});
            }
            if (boundary.y) |y| {
                std.debug.print("[INFO] y = {d}\n", .{y});
            }
            if (boundary.z) |z| {
                std.debug.print("[INFO] z = {d}\n", .{z});
            }
            std.debug.print("[INFO]\n", .{});
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

    fn display(self: *const Self) void {
        std.debug.print("[INFO] CRD\n", .{});
        for (self.coordinates) |coordinate| {
            std.debug.print("[INFO] id = {d}, x = {d}, y = {d}, z = {d}\n", .{ coordinate.id, coordinate.x, coordinate.y, coordinate.z });
        }
        std.debug.print("[INFO]\n", .{});
    }
};

const Class = struct {
    name: []const u8,
    mass: f64,
    charge: f64,
};

const LJ = struct {
    class1: []const u8,
    class2: []const u8,
    epsilon: f64,
    sigma: f64,
};

const Par = struct {
    classes: []const Class,
    lj: []const LJ,

    const Self = @This();

    fn display(self: *const Self) void {
        std.debug.print("[INFO] CLASSES\n", .{});
        for (self.classes) |class| {
            std.debug.print("[INFO] name = {s}, mass = {d}, charge = {d}\n", .{ class.name, class.mass, class.charge });
        }
        std.debug.print("[INFO]\n", .{});

        std.debug.print("[INFO] LJ\n", .{});
        for (self.lj) |lj| {
            std.debug.print("[INFO] class1 = {s}, class2 = {s}, epsilon = {d}, sigma = {d}\n", .{ lj.class1, lj.class2, lj.epsilon, lj.sigma });
        }
        std.debug.print("[INFO]\n", .{});
    }
};

const Particle = struct {
    id: u64,
    class: []const u8,
};

const Top = struct {
    particles: []Particle,

    const Self = @This();

    fn display(self: *const Self) void {
        std.debug.print("[INFO] TOP\n", .{});
        for (self.particles) |particle| {
            std.debug.print("[INFO] id = {d}, class = {s}\n", .{ particle.id, particle.class });
        }
        std.debug.print("[INFO]\n", .{});
    }
};

pub fn main() !void {
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
    cfg.display();

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
    crd.display();

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
    par.display();

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
    top.display();
}
