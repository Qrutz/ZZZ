package bouncing_balls;

/**
 * The physics model.
 * 
 * This class is where you should implement your bouncing balls model.
 * 
 * The code has intentionally been kept as simple as possible, but if you wish,
 * you can improve the design.
 * 
 * @author Simon Robillard
 *
 */
class Model {

	double areaWidth, areaHeight;

	Ball[] balls;

	// gravity constant
	private static final double GRAVITY = 9.81;

	// how much energy gets lost when hitting walls
	private static final double DAMPING = 0.95;

	Model(double width, double height) {
		areaWidth = width;
		areaHeight = height;

		// Initialize the model with a few balls
		balls = new Ball[2];
		balls[0] = new Ball(width / 3, height * 0.9, 1.2, 1.6, 0.2);
		balls[1] = new Ball(2 * width / 3, height * 0.7, -0.6, 0.6, 0.3);
	}

	void step(double deltaT) {
		// move balls forward in time
		for (Ball b : balls) {
			// gravity pulls balls down
			b.vy += deltaT * (-GRAVITY);

			// update position based on current speed
			b.x += deltaT * b.vx;
			b.y += deltaT * b.vy;
		}

		// check if balls are hitting each other
		for (int i = 0; i < balls.length; i++) {
			for (int j = i + 1; j < balls.length; j++) {
				if (ballsColliding(balls[i], balls[j])) {
					handleBallCollision(balls[i], balls[j]);
				}
			}
		}

		// check if balls hit the walls
		for (Ball b : balls) {
			// hit left or right wall
			if (b.x < b.radius || b.x > areaWidth - b.radius) {
				b.vx *= -DAMPING;
				// keep ball inside the area
				if (b.x < b.radius)
					b.x = b.radius;
				if (b.x > areaWidth - b.radius)
					b.x = areaWidth - b.radius;
			}
			// hit top or bottom wall
			if (b.y < b.radius || b.y > areaHeight - b.radius) {
				b.vy *= -DAMPING;
				// keep ball inside the area
				if (b.y < b.radius)
					b.y = b.radius;
				if (b.y > areaHeight - b.radius)
					b.y = areaHeight - b.radius;
			}
		}
	}

	// check if two balls are touching
	private boolean ballsColliding(Ball ball1, Ball ball2) {
		double dx = ball2.x - ball1.x;
		double dy = ball2.y - ball1.y;
		double distance = Math.sqrt(dx * dx + dy * dy);
		return distance <= (ball1.radius + ball2.radius);
	}

	// make balls bounce off each other
	private void handleBallCollision(Ball ball1, Ball ball2) {
		// figure out which direction the collision is happening
		double dx = ball2.x - ball1.x;
		double dy = ball2.y - ball1.y;
		double distance = Math.sqrt(dx * dx + dy * dy);

		// get the collision angle using polar coordinates
		double[] polar_pos = rectToPolar(dx, dy);
		double collision_angle = polar_pos[1];

		// convert velocities to polar coordinates
		double[] polar_v1 = rectToPolar(ball1.vx, ball1.vy);
		double[] polar_v2 = rectToPolar(ball2.vx, ball2.vy);

		// get velocity components along the collision line
		double v1_collision = polar_v1[0] * Math.cos(polar_v1[1] - collision_angle);
		double v2_collision = polar_v2[0] * Math.cos(polar_v2[1] - collision_angle);

		// if they're moving apart, don't do anything
		if (v2_collision - v1_collision > 0)
			return;

		// apply collision physics - relative velocity reverses sign
		double relative_vel = -(v2_collision - v1_collision);

		// conservation of momentum
		double total_momentum = ball1.mass * v1_collision + ball2.mass * v2_collision;

		// solve for new velocities
		double v1_new = (total_momentum - ball2.mass * relative_vel) / (ball1.mass + ball2.mass);
		double v2_new = v1_new + relative_vel;

		// keep the tangential velocity components (perpendicular to collision)
		double[] tan_v1 = polarToRect(polar_v1[0] * Math.sin(polar_v1[1] - collision_angle),
				collision_angle + Math.PI / 2);
		double[] tan_v2 = polarToRect(polar_v2[0] * Math.sin(polar_v2[1] - collision_angle),
				collision_angle + Math.PI / 2);

		// combine tangential and collision components
		ball1.vx = tan_v1[0] + v1_new * Math.cos(collision_angle);
		ball1.vy = tan_v1[1] + v1_new * Math.sin(collision_angle);

		ball2.vx = tan_v2[0] + v2_new * Math.cos(collision_angle);
		ball2.vy = tan_v2[1] + v2_new * Math.sin(collision_angle);

		// lose a tiny bit of energy when balls hit each other
		double damping_factor = 0.995;
		ball1.vx *= damping_factor;
		ball1.vy *= damping_factor;
		ball2.vx *= damping_factor;
		ball2.vy *= damping_factor;

		// push balls apart so they don't get stuck
		double overlap = (ball1.radius + ball2.radius) - distance;
		double[] separation_polar = polarToRect(overlap / 2, collision_angle);

		ball1.x -= separation_polar[0];
		ball1.y -= separation_polar[1];
		ball2.x += separation_polar[0];
		ball2.y += separation_polar[1];
	}

	// convert rectangular coordinates to polar
	private double[] rectToPolar(double x, double y) {
		double r = Math.sqrt(x * x + y * y);
		double theta = Math.atan2(y, x);
		return new double[] { r, theta };
	}

	// convert polar coordinates to rectangular
	private double[] polarToRect(double r, double theta) {
		double x = r * Math.cos(theta);
		double y = r * Math.sin(theta);
		return new double[] { x, y };
	}

	// class to represent each ball
	class Ball {

		Ball(double x, double y, double vx, double vy, double r) {
			this.x = x;
			this.y = y;
			this.vx = vx;
			this.vy = vy;
			this.radius = r;
			this.mass = r; // mass proportional to radius
		}

		// ball properties
		double x, y, vx, vy, radius, mass;
	}
}
