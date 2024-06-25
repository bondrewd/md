const std = @import("std");
const F = @import("precision.zig").F;

pub const Class = struct {
    name: []const u8,
    mass: F,
    charge: F,
};

pub const LJ = struct {
    name: []const u8,
    epsilon: F,
    sigma: F,
};

pub const Par = struct {
    classes: []const Class,
    lj: []const LJ,

    const Self = @This();

    pub fn writeLog(self: *const Self, writer: std.fs.File.Writer) !void {
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
