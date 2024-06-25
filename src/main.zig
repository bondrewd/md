const std = @import("std");
const builtin = @import("builtin");
const F = @import("precision.zig").F;
const BOLTZMANN: F = @import("constants.zig").BOLTZMANN;
const Cfg = @import("cfg.zig").Cfg;
const System = @import("system.zig").System;
const Crd = @import("crd.zig").Crd;
const Par = @import("par.zig").Par;
const Top = @import("top.zig").Top;
const XYZFile = @import("xyz.zig").XYZFile;
const CSVFile = @import("csv.zig").CSVFile;

pub fn main() !void {
    const stderr = std.io.getStdErr().writer();

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        const deinit_status = gpa.deinit();
        if (deinit_status == .leak) std.testing.expect(false) catch @panic("TEST FAIL");
    }
    const allocator = gpa.allocator();

    var env_map = try std.process.getEnvMap(allocator);
    defer env_map.deinit();

    std.debug.print("[INFO] MD\n", .{});
    std.debug.print("[INFO] cpu = {s}\n", .{@tagName(builtin.target.cpu.arch)});
    std.debug.print("[INFO] os = {s}\n", .{@tagName(builtin.target.os.tag)});
    std.debug.print("[INFO] user = {s}\n", .{env_map.get("USER") orelse "unknown"});
    std.debug.print("[INFO] zig = {s}\n", .{builtin.zig_version_string});
    std.debug.print("[INFO] mode = {s}\n", .{@tagName(builtin.mode)});
    std.debug.print("[INFO] precision = {s}\n", .{@typeName(F)});
    std.debug.print("[INFO]\n", .{});

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
    try cfg.writeLog(stderr);

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
    try crd.writeLog(stderr);

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
    try par.writeLog(stderr);

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
    try top.writeLog(stderr);

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
        // Open XYZ file
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

        // Open CSV file
        var csv_file_path: ?[]u8 = null;
        if (cfg.output) |output| {
            if (output.csv) |csv| {
                if (std.fs.path.dirname(cfg_file_path)) |dir| {
                    csv_file_path = try std.fs.path.join(allocator, &[_][]const u8{ dir, csv });
                } else {
                    csv_file_path = try std.fs.path.join(allocator, &[_][]const u8{ ".", csv });
                }
            }
        }
        defer if (csv_file_path) |path| allocator.free(path);
        const csv_file = CSVFile.init(csv_file_path);
        defer csv_file.deinit();

        // Initial output
        try csv_file.writeHeader();
        try csv_file.writeState(&system, 0);
        try xyz_file.writeFrame(&system);

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

            if (step % 100 == 0) try csv_file.writeState(&system, step);
            if (step % 100 == 0) try xyz_file.writeFrame(&system);
        }
    }
}
