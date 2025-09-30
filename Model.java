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

	// gravity pulls the balls down
	private static final double GRAVITY = 9.81;

	// how much energy is lost when hitting walls (0.9 means lose 10% each bounce)
	private static final double DAMPING = 0.9;

	Model(double width, double height) {
		areaWidth = width;
		areaHeight = height;

		// Initialize the model with a few balls
		balls = new Ball[2];
		balls[0] = new Ball(width / 3, height * 0.9, 1.2, 1.6, 0.2);
		balls[1] = new Ball(2 * width / 3, height * 0.7, -0.6, 0.6, 0.3);
	}

	void step(double deltaT) {
		// move the balls forward in time
		for (Ball b : balls) {
			// gravity makes them fall faster
			b.vy += deltaT * (-GRAVITY);

			// update position based on current speed
			b.x += deltaT * b.vx;
			b.y += deltaT * b.vy;
		}

		// check if any balls are hitting each other
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
				b.vx *= -DAMPING; // bounce back and lose some energy
				// make sure they don't get stuck in the wall
				if (b.x < b.radius)
					b.x = b.radius;
				if (b.x > areaWidth - b.radius)
					b.x = areaWidth - b.radius;
			}
			// hit top or bottom wall
			if (b.y < b.radius || b.y > areaHeight - b.radius) {
				b.vy *= -DAMPING; // bounce back and lose some energy
				// make sure they don't get stuck in the wall
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

	// make balls bounce off each other realistically
	private void handleBallCollision(Ball ball1, Ball ball2) {
		// figure out which direction the collision is happening
		double dx = ball2.x - ball1.x;
		double dy = ball2.y - ball1.y;
		double distance = Math.sqrt(dx * dx + dy * dy);

		// normalize the collision direction
		double nx = dx / distance;
		double ny = dy / distance;

		// check if balls are actually moving towards each other
		double dvx = ball2.vx - ball1.vx;
		double dvy = ball2.vy - ball1.vy;
		double dvn = dvx * nx + dvy * ny;

		// if they're moving apart, don't do anything
		if (dvn > 0)
			return;

		// use the physics from collision.txt to calculate new velocities
		// this makes the collision look realistic

		// remember the old velocities
		double u1x = ball1.vx;
		double u1y = ball1.vy;
		double u2x = ball2.vx;
		double u2y = ball2.vy;

		// calculate how fast they were approaching each other
		double R = dvn;

		// calculate total momentum in the collision direction
		double I = ball1.mass * u1x * nx + ball1.mass * u1y * ny +
				ball2.mass * u2x * nx + ball2.mass * u2y * ny;

		// solve the physics equations to get new velocities
		double v1_new = (I - ball2.mass * (-R)) / (ball1.mass + ball2.mass);
		double v2_new = v1_new + (-R);

		// update the ball velocities
		ball1.vx = u1x + (v1_new - (u1x * nx + u1y * ny)) * nx;
		ball1.vy = u1y + (v1_new - (u1x * nx + u1y * ny)) * ny;

		ball2.vx = u2x + (v2_new - (u2x * nx + u2y * ny)) * nx;
		ball2.vy = u2y + (v2_new - (u2x * nx + u2y * ny)) * ny;

		// lose a tiny bit of energy when balls hit each other
		double damping_factor = 0.99;
		ball1.vx *= damping_factor;
		ball1.vy *= damping_factor;
		ball2.vx *= damping_factor;
		ball2.vy *= damping_factor;

		// push balls apart so they don't get stuck
		double overlap = (ball1.radius + ball2.radius) - distance;
		double separationX = overlap * nx / 2;
		double separationY = overlap * ny / 2;

		ball1.x -= separationX;
		ball1.y -= separationY;
		ball2.x += separationX;
		ball2.y += separationY;
	}

	// class to represent each ball
	class Ball {

		Ball(double x, double y, double vx, double vy, double r) {
			this.x = x;
			this.y = y;
			this.vx = vx;
			this.vy = vy;
			this.radius = r;
			this.mass = Math.PI * r * r; // bigger balls are heavier
		}

		// ball properties - position, velocity, size, and weight
		double x, y, vx, vy, radius, mass;
	}
}
