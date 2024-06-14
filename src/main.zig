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

const Configuration = struct {
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
    x: f64,
    y: f64,
    z: f64,
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

    const parsed = try std.json.parseFromSlice(Configuration, allocator, cfg_data, .{});
    defer parsed.deinit();

    const cfg = parsed.value;
    cfg.display();

    const cfg_file_dir = std.fs.path.dirname(cfg_file_path);
    var crd_file_path: []u8 = undefined;
    if (cfg_file_dir) |dir| {
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

    var crd = std.AutoHashMap(u64, Coordinate).init(allocator);
    defer crd.deinit();

    var lines = std.mem.splitSequence(u8, crd_data, "\n");
    while (lines.next()) |line| {
        if (line.len == 0) continue;
        if (std.mem.startsWith(u8, line, "#")) continue;
        var tokens = std.mem.splitSequence(u8, line, ",");
        const i = try std.fmt.parseInt(u64, tokens.next().?, 10);
        const x = try std.fmt.parseFloat(f64, tokens.next().?);
        const y = try std.fmt.parseFloat(f64, tokens.next().?);
        const z = try std.fmt.parseFloat(f64, tokens.next().?);
        const old = try crd.fetchPut(i, Coordinate{ .x = x, .y = y, .z = z });
        if (old) |_| {
            std.debug.print("[ERROR] Duplicate index {d} in crd file\n", .{i});
            return;
        }
    }

    var entries = crd.iterator();
    while (entries.next()) |entry| {
        const i = entry.key_ptr.*;
        const x = entry.value_ptr.x;
        const y = entry.value_ptr.y;
        const z = entry.value_ptr.z;
        std.debug.print("[INFO] i = {d}, x = {d}, y = {d}, z = {d}\n", .{ i, x, y, z });
    }
}
