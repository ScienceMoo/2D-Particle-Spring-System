package comp559.particle;

import java.io.Serializable;
import java.util.ArrayList;

import javax.vecmath.Color3f;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * Particle class that contains particle properties (e.g., mass), 
 * initial positions and velocities, current position and velocities 
 * and a force accumulator for computing the total force acting on this particle.
 * @author kry
 */
public class Particle implements Serializable {
    
    /** Identifies this particles position in the particle list */
    int index;
    
    public boolean pinned = false;

    public boolean visible = true;
    
    Color3f color = new Color3f(0,.95f,0);
    
    float size = 10;
    
    double mass = 1;

    /** Position of the particle */
    Point3d p = new Point3d();

    /** Velocity of the particle */
    Vector3d v = new Vector3d();

    /** Initial position of the particle */
    Point3d p0 = new Point3d();

    /** Initial velocity of the particle */
    Vector3d v0 = new Vector3d();

    /**
     * Forces acting on the particle.  Accumulate them here, but
     * remember to reset them to zero before starting.
     */
    Vector3d f = new Vector3d();

    /**
     * A list of springs that use this particle.  This list is only needed
     * to adjust rest lengths when dragging particles around.
     * This is only used for UI... it is probably not needed for simulation 
     */
    ArrayList<Spring> springs = new ArrayList<Spring>();
    
    /**
     * Creates a particle with the given position and velocity
     * @param p0
     * @param v0
     */
    public Particle( Point3d p0, Vector3d v0 ) {
        this.p0.set(p0);
        this.v0.set(v0);
        reset();
    }

    /**
     * Creates a particle with the given position and velocity and mass
     * @param p0
     * @param v0
     * @param p_mass
     */
    public Particle( Point3d p0, Vector3d v0, double p_mass) {
        this.p0.set(p0);
        this.v0.set(v0);
        //TODO: dpi stuff
//        mass = p_mass * 2;
        mass = p_mass;
        reset();
    }

    /**
     * Constructor for invisible particles
     * @param p0
     * @param v0
     * @param p_mass
     * @param vis
     */
    public Particle( Point3d p0, Vector3d v0, double p_mass, boolean vis) {
        this.p0.set(p0);
        this.v0.set(v0);
        //TODO: dpi stuff
        mass = p_mass;
//        mass = p_mass * 2;
        visible = vis;
        pinned = true;
        reset();
    }
    
    /**
     * Resets the position of this particle
     */
    public void reset() {
        p.set(p0);
        v.set(v0);
        f.set(0,0,0);
    }
    
    /**
     * Clears all forces acting on this particle
     */
    public void clearForce() {

        f.set(0,0,0);
    }
    
    /**
     * Adds the given force to this particle
     * @param force
     */
    public void addForce( Vector3d force ) {

        f.add(force);
    }
    
    /**
     * Computes the distance of a point to this particle
     * @param x
     * @param y
     * @return the distance
     */
    public double distance( double x, double y, double z ) {
        Point3d tmp = new Point3d( x, y, z );
        return tmp.distance(p);
    }
    
    @Override
    public String toString() {
        return f.toString();
    }
    
}
