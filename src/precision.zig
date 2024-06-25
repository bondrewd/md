const config = @import("config");

pub const F = switch (config.precision) {
    .Single => f32,
    .Double => f64,
};
