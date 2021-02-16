package comp559.particle;

import javax.vecmath.Vector2d;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

import java.io.Serializable;

/**
 * Spring class for 599 assignment 1
 * @author kry
 */
public class Spring implements Serializable {

    Particle p1 = null;
    Particle p2 = null;
    
    /** Spring stiffness, sometimes written k_s in equations */
    public static double k = 1;
    /** Spring damping (along spring direction), sometimes written k_d in equations */
    public static double c = 1;
    /** Rest length of this spring */
    double l0 = 0;
    
    /**
     * Creates a spring between two particles
     * @param p1
     * @param p2
     */
    public Spring( Particle p1, Particle p2 ) {
        this.p1 = p1;
        this.p2 = p2;
        recomputeRestLength();
        p1.springs.add(this);
        p2.springs.add(this);
        this.apply();
    }

    /**
     * Computes and sets the rest length based on the original position of the two particles 
     */
    public void recomputeRestLength() {

        l0 = p1.p0.distance( p2.p0 );
    }
    
    /**
     * Applies the spring force by adding a force to each particle
     */
    public void apply() {
        // TODO: Objective 1, FINISH THIS CODE!
        Vector2d force = new Vector2d();
//        System.out.print("force");
//        System.out.println(force);
//        System.out.print("p2.p");
//        System.out.println(p2.p);
//        System.out.print("p1.p");
//        System.out.println(p1.p);

        // calculate force from spring deformation
        force.sub( p2.p, p1.p );
//        System.out.print("force after sub");
//        System.out.println(force);

        double l = force.length();
//        System.out.print("l");
//        System.out.println(l);

        force.normalize();
//        System.out.print("force after normalize");
//        System.out.println(force);

        // stiffness determines how hard the spring will try to go back to resting length
        force.scale( (l-l0)*k );
        p1.addForce(force);
        force.scale(-1);
        p2.addForce(force);

        // calculate force from kinetic energy
        force.sub( p2.p, p1.p );
        force.normalize();
        Vector2d v = new Vector2d();
        v.sub(p2.v, p1.v);
        double rv = force.dot(v);

        // damping determines how much kinetic energy will be lost due to friction
        force.scale(c*rv);
        p1.addForce(force);
        force.scale(-1);
        p2.addForce(force);
    }
   
    /** the functions below are for the backwards Euler solver */
    
    /**
     * Computes the force and adds it to the appropriate components of the force vector.
     * (This function is something you might use for a backward Euler integrator)
     * @param f
     */
    public void addForce( Vector f ) {
        // TODO: Objective 8, FINISH THIS CODE for backward Euler method (probably very simlar to what you did above)
        Vector2d force = new Vector2d();

        force.sub( p2.p, p1.p );
        double l = force.length();
        force.normalize();
        force.scale( (l-l0)*k );
        f.add(p1.index*3,force.x);
        f.add(p1.index*3+1,force.y);
        force.scale(-1);
        f.add(p2.index*3,force.x);
        f.add(p2.index*3+1,force.y);

        // calculate force from kinetic energy
        force.sub( p2.p, p1.p );
        force.normalize();
        Vector2d v = new Vector2d();
        v.sub(p2.v, p1.v);
        double rv = force.dot(v);
        force.scale(c*rv);
        f.add(p1.index*3,force.x);
        f.add(p1.index*3+1,force.y);
        force.scale(-1);
        f.add(p2.index*3,force.x);
        f.add(p2.index*3+1,force.y);
    }

    public void addParticle(Matrix mtrx, Particle p1, Particle p2, DenseMatrix tmpM) {
        mtrx.add(p1.index*3, p2.index*3, tmpM.get(0, 0));
        mtrx.add(p1.index*3, p2.index*3+1, tmpM.get(0, 1));
        mtrx.add(p1.index*3+1, p2.index*3, tmpM.get(1, 0));
        mtrx.add(p1.index*3+1, p2.index*3+1, tmpM.get(1, 1));
    }
    
    /**
     * Adds this springs contribution to the stiffness matrix
     * @param dfdx
     */
    public void addDfdx( Matrix dfdx ) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward euler integration
        // dfdx = -k * ((1 - l0/|l|)(I - l * l^T) + l * l^T)

        DenseVector tmpV = new DenseVector(2);
        DenseMatrix tmpM = new DenseMatrix(2,2);
        DenseMatrix I = new DenseMatrix(2,2);

        Vector2d l = new Vector2d(p1.p.x - p2.p.x, p1.p.y - p2.p.y);

        tmpV.set(0,l.x/l.length());
        tmpV.set(1,l.y/l.length());

        tmpM.rank1(tmpV);
        I.set(tmpM);
        I.scale(-1);
        I.add(0, 0, 1.0);
        I.add(1, 1, 1.0);

        I.scale(1 - l0/l.length());

        tmpM.add(I);
        tmpM.scale(-k);

        addParticle(dfdx, p1, p1, tmpM);
        addParticle(dfdx, p2, p2, tmpM);

        tmpM.scale(-1);

        addParticle(dfdx, p1, p2, tmpM);
        addParticle(dfdx, p2, p1, tmpM);

    }   
 
    /**
     * Adds this springs damping contribution to the implicit damping matrix
     * @param dfdv
     */
    public void addDfdv( Matrix dfdv ) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward Euler integration
        // dfdv = -c * l * l^T

        DenseVector tmpV = new DenseVector(2);
        DenseMatrix tmpM = new DenseMatrix(2,2);

        Vector2d l = new Vector2d(p1.p.x - p2.p.x, p1.p.y - p2.p.y);

        tmpV.set(0,l.x/l.length());
        tmpV.set(1,l.y/l.length());

        tmpM.rank1(tmpV);
        tmpM.scale(-c);

        addParticle(dfdv, p1, p1, tmpM);
        addParticle(dfdv, p2, p2, tmpM);

        tmpM.scale(-1);

        addParticle(dfdv, p1, p2, tmpM);
        addParticle(dfdv, p2, p1, tmpM);
    } 
    
}
