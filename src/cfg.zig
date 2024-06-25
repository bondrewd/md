const std = @import("std");
const builtin = @import("builtin");

const F = @import("precision.zig").F;
const BOLTZMANN: F = 1.987204259e-3;

const InputBlock = struct {
    crd: []u8,
    par: []u8,
    top: []u8,

    const Self = @This();

    fn writeLog(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] INPUT\n", .{});
        try writer.print("[INFO] crd = {s}\n", .{self.crd});
        try writer.print("[INFO] par = {s}\n", .{self.par});
        try writer.print("[INFO] top = {s}\n", .{self.top});
    }
};

const OutputBlock = struct {
    crd: ?[]u8 = null,
    xyz: ?[]u8 = null,
    log: ?[]u8 = null,
    csv: ?[]u8 = null,

    const Self = @This();

    fn writeLog(self: *const Self, writer: std.fs.File.Writer) !void {
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
        if (self.log) |log| {
            try writer.print("[INFO] log = {s}\n", .{log});
        } else {
            try writer.print("[INFO] log =\n", .{});
        }
        if (self.csv) |csv| {
            try writer.print("[INFO] csv = {s}\n", .{csv});
        } else {
            try writer.print("[INFO] csv =\n", .{});
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
    dt: F,
    steps: u64,
    ensemble: ?EnsembleKind = null,
    temperature: ?F = null,
    pressure: ?F = null,
    thermostat: ?ThermostatKind = null,
    barostat: ?BarostatKind = null,

    const Self = @This();

    fn writeLog(self: *const Self, writer: std.fs.File.Writer) !void {
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

    fn writeLog(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] RNG\n", .{});
        try writer.print("[INFO] seed = {d}\n", .{self.seed});
    }
};

const BoundaryBlock = struct {
    x: ?F = null,
    y: ?F = null,
    z: ?F = null,

    const Self = @This();

    fn writeLog(self: *const Self, writer: std.fs.File.Writer) !void {
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

pub const Cfg = struct {
    input: InputBlock,
    output: ?OutputBlock = null,
    dynamics: ?DynamicsBlock = null,
    boundary: ?BoundaryBlock = null,
    rng: ?Rng = null,

    const Self = @This();

    pub fn writeLog(self: *const Self, writer: std.fs.File.Writer) !void {
        try self.input.writeLog(writer);
        try writer.print("[INFO] \n", .{});
        if (self.output) |output| {
            try output.writeLog(writer);
            try writer.print("[INFO] \n", .{});
        }
        if (self.dynamics) |dynamics| {
            try dynamics.writeLog(writer);
            try writer.print("[INFO] \n", .{});
        }
        if (self.rng) |rng| {
            try rng.writeLog(writer);
            try writer.print("[INFO] \n", .{});
        }
        if (self.boundary) |boundary| {
            try boundary.writeLog(writer);
            try writer.print("[INFO] \n", .{});
        }
    }
};
