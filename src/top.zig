const std = @import("std");
const Particle = @import("particle.zig").Particle;

pub const Top = struct {
    particles: []Particle,

    const Self = @This();

    pub fn writeLog(self: *const Self, writer: std.fs.File.Writer) !void {
        try writer.print("[INFO] TOP\n", .{});
        try writer.print("[INFO] Number of particles: {d}\n", .{self.particles.len});
        try writer.print("[INFO]\n", .{});
    }
};
