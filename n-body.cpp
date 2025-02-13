#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>      // for output saving.
#include <cmath>        // for sqrt.
#include <chrono>       // for timing.

const double G = 6.674e-11;
const double SOFTENING = 0.0001;

// Planets Masses
const double SUN_MASS = 1.989e30;
const double MERCURY_MASS = 3.285e23;
const double VENUS_MASS = 4.867e24;
const double EARTH_MASS = 5.972e24;
const double MARS_MASS = 6.39e23;
const double JUPITER_MASS = 1.898e27;
const double SATURN_MASS = 5.683e26;
const double URANUS_MASS = 8.681e25;
const double NEPTUNE_MASS = 1.024e26;
const double PLANET_MASSES[] = {SUN_MASS, MERCURY_MASS, VENUS_MASS, EARTH_MASS, MARS_MASS, JUPITER_MASS, SATURN_MASS, URANUS_MASS, NEPTUNE_MASS};

struct Particle {
    double mass;            // mass
    double x, y, z;         // positions
    double vx, vy, vz;      // velocities
    double fx, fy, fz;      // forces
};

void initialize_particles(std::vector<Particle>& particles, int num_particles) {
    for (size_t i = 0; i < num_particles; i++) {
        // use a for loop to initialize all particles with planets masses & position at 0 velocity & force.
        particles.push_back({
            PLANET_MASSES[std::rand() % 9],                      // random planet mass
            rand() % 10000,              // random position x
            rand() % 10000,              // random position y.
            rand() % 10000,              // random position z.
            0.0, 0.0, 0.0,                                  // velocities
            0.0, 0.0, 0.0                                   // forces.
            });
        }
    }    

void compute_forces(std::vector<Particle>& particles) {
    // Reset the Forces to 0.
    for (size_t i = 0; i < particles.size(); i++) {
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        particles[i].fz = 0.0;
    }

    // 2 Loops for to calculate every pair.
    for (size_t i = 0; i < particles.size(); i++) {
        for (size_t j = i+1; j < particles.size(); j++) {
            // Distance: Particle 1 - Particle 2
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double dz = particles[j].z - particles[i].z;

            // r^2 = Distance Squared
            double distance_squared = dx * dx + dy * dy + dz * dz + SOFTENING;
            double dist = std::sqrt(distance_squared); // mangnitude of distance aka || l1 - l2 ||

            // Force = Gravity ((Mass 1 * Mass 2) / r^2) aka Newton's Law
            double force = (G * particles[i].mass * particles[j].mass) / distance_squared;

            // Force of X/Y/Z = F * Distance Vector / distance (r^2) along the X/Y/Z axis
            double fx = force * dx / dist;
            double fy = force * dy / dist;
            double fz = force * dz / dist;

            // The direction of the force is from the particle that applies the force and towards the particle that it is applied on.
            // Newton's Third Law
            particles[i].fx += fx;
            particles[i].fy += fy;
            particles[i].fz += fz;

            particles[j].fx -= fx;
            particles[j].fy -= fy;
            particles[j].fz -= fz;
        }
    }
}

void update_particles(std::vector<Particle>& particles, double dt) {
    for (auto& p : particles) {
        // Acceleration = Force / Mass.
        double ax = p.fx / p.mass;
        double ay = p.fy / p.mass;
        double az = p.fz / p.mass;

        // Velocity = Old Velocity + Acceleration * Time.
        p.vx += ax * dt;
        p.vy += ay * dt;
        p.vz += az * dt;

        // Position = Old Position + Velocity * Time.
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
    }
}

void save_state(const std::vector<Particle>& particles, std::ofstream& file) {
    file << particles.size();
    for (const auto& p : particles) {
        file << "\t" << p.mass << "\t" << p.x << "\t" << p.y << "\t" << p.z;
        file << "\t" << p.vx << "\t" << p.vy << "\t" << p.vz;
        file << "\t" << p.fx << "\t" << p.fy << "\t" << p.fz;
    }
    file << "\n";
}


int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << "<num_particles> <dt> <iterations> <output_interval>\n";
        return 1;
    }

    int num_particles = std::stoi(argv[1]);
    double dt = std::stod(argv[2]);
    int iterations = std::stoi(argv[3]);
    int output_interval = std::stoi(argv[4]);

    std::vector<Particle> particles;     // declare particles.
    initialize_particles(particles, num_particles);     // initalization.
    std::ofstream file("output.tsv");                        // create an output file.

    auto start = std::chrono::steady_clock::now();    // start time.

    for (size_t step = 0; step < iterations; step++) {
        compute_forces(particles);
        update_particles(particles, dt);   

        if (step % output_interval == 0) {      // Save the output by interval.
            save_state(particles, file);
        }
    }
    auto finish = std::chrono::steady_clock::now();      // end time.
    std::chrono::duration<double> elapsed_seconds = finish - start;
    std::cout << "Simulation Time: " << elapsed_seconds.count() << "seconds.\n";

    file.close();   // close file.
    return 0;
}