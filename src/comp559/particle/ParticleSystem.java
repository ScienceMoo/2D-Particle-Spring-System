package comp559.particle;

import java.util.LinkedList;
import java.util.List;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.SceneGraphNode;

import java.io.*;

/**
 * Implementation of a simple particle system
 * @author kry
 */
public class ParticleSystem implements SceneGraphNode, Function, Filter {

    public class ParticleSystemInfo implements Serializable {
        public List<Particle> particles;
        public List<Spring> springs;

        public ParticleSystemInfo() {
            particles = new LinkedList<Particle>();
            springs = new LinkedList<Spring>();
        }
    }

    private ParticleSystemInfo info;
    
    /**
     * Creates an empty particle system
     */
    public ParticleSystem() {
        info = new ParticleSystemInfo();
    }

    /**
     * Create a saved test system from the file
     * @param filename
     */
    public void loadSystem( String filename ) {
        // Deserialization
        try {

            // Reading the object from a file
            FileInputStream file = new FileInputStream(filename);
            ObjectInputStream in = new ObjectInputStream(file);

            // Method for deserialization of object
            ParticleSystemInfo loadedSystem = (ParticleSystemInfo) in.readObject();
            this.info.particles = loadedSystem.particles;
            this.info.springs = loadedSystem.springs;

            in.close();
            file.close();
        }

        catch (IOException ex) {
            System.out.println("IOException is caught");
        }

        catch (ClassNotFoundException ex) {
            System.out.println("ClassNotFoundException" +
                    " is caught");
        }
    }

//    private void writeObject(ObjectOutputStream out) throws IOException {
//        out.defaultWriteObject();
//    }
//    private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
//        in.defaultReadObject();
//    }

    /**
     * Saves current particle system to a file that can be loaded later
     * @param filename
     */
    public void saveSystem( String filename ) {
        // Serialization
        try {

            // create new file
            File myObj = new File("/Users/moo/comp559W19/savedSystems/" + filename);

            // Saving of object in a file
            FileOutputStream file = new FileOutputStream("/Users/moo/comp559W19/savedSystems/" + filename);
            ObjectOutputStream out = new ObjectOutputStream(file);

            // Method for serialization of object
            out.writeObject(this.info);
            out.flush();
            out.close();
//            file.close();

        }

        catch (IOException ex) {
            System.out.println("IOException is caught");
            System.out.println(ex.getCause());
            System.out.println(ex.getStackTrace());
            System.out.println(ex.getMessage());
            System.out.println("IOException is caught 2");
        }
    }

    /**
     * Creates one of a number of simple test systems.
     * @param which
     */
    public void createSystem( int which ) {
        
        if ( which == 1) {        
            Point2d p = new Point2d( 100, 100 );
            Vector2d d = new Vector2d( 20, 0 );            
            Particle p1, p2, p3, p4;
            p1 = new Particle( p.x - d.y, p.y + d.x, 0, 0 );
            p1.index = info.particles.size();
            info.particles.add( p1 );
            p2 = new Particle( p.x + d.y, p.y - d.x, 0, 0 );
            p2.index = info.particles.size();
            info.particles.add( p2 );
            info.springs.add( new Spring ( p1, p2 ) );
            p1.pinned = true;
            p2.pinned = true;            
            p.add( d );
            p.add( d );                    
            int N = 10;
            for (int i = 1; i < N; i++ ) {                
                //d.set( 20*Math.cos(i*Math.PI/N), 20*Math.sin(i*Math.PI/N) );                
                p3 = new Particle( p.x - d.y, p.y + d.x, 0, 0 );
                p3.index = info.particles.size();
                info.particles.add( p3 );
                p4 = new Particle( p.x + d.y, p.y - d.x, 0, 0 );
                p4.index = info.particles.size();
                info.particles.add( p4 );
                info.springs.add( new Spring ( p3, p1 ) );
                info.springs.add( new Spring ( p3, p2 ) );
                info.springs.add( new Spring ( p4, p1 ) );
                info.springs.add( new Spring ( p4, p2 ) );
                info.springs.add( new Spring ( p4, p3 ) );
                p1 = p3;
                p2 = p4;                
                p.add( d );
                p.add( d );            
            }
        } else if ( which == 2) {
            Particle p1 = new Particle( 320, 100, 0, 0 );
            p1.index = info.particles.size();
            info.particles.add( p1 );
            Particle p2 = new Particle( 320, 200, 0, 0 );
            p2.index = info.particles.size();
            info.particles.add( p2 );
            p1.pinned = true;
            info.springs.add( new Spring( p1, p2 ) );
        } else if ( which == 3 ) {
            int ypos = 100;
            Particle p0 = null;
            Particle p1, p2;
            p1 = new Particle( 320, ypos, 0, 0 );
            p1.index = info.particles.size();
            p1.pinned = true;
            info.particles.add( p1 );
            int N = 10;
            for ( int i = 0; i < N; i++ ) {
                ypos += 20;
                p2 = new Particle( 320, ypos, 0, 0 );
                p2.index = info.particles.size();
                info.particles.add( p2 );
                info.springs.add( new Spring( p1, p2 ) );
                // Hum.. this is not great in comparison to a proper bending energy...
                // use Maple to generate some code though, as it is painful to write by hand! :(
                if ( p0 != null ) info.springs.add( new Spring( p2, p0 ) );
                p0 = p1;
                
                p1 = p2;
            }
        }
    }
    
    /**
     * Gets the particles in the system
     * @return the particle set
     */
    public List<Particle> getParticles() {
        return info.particles;
    }
    
    /**
     * Gets the springs in the system
     * @return the spring list
     */
    public List<Spring> getSprings() {
    	return info.springs;
    }
    
    /**
     * Resets the positions of all particles
     */
    public void resetParticles() {
        for ( Particle p : info.particles ) {
            p.reset();
        }
        time = 0;
    }
    
    /**
     * Deletes all particles
     */
    public void clearParticles() {
        info.particles.clear();
        info.springs.clear();
    }
    
    /**
     * Gets the phase space state of the particle system
     * @param phaseSpaceState
     */
    public void getPhaseSpace( double[] phaseSpaceState ) {
        int count = 0;
        for ( Particle p : info.particles ) {
            phaseSpaceState[count++] = p.p.x;
            phaseSpaceState[count++] = p.p.y;
            phaseSpaceState[count++] = p.v.x;
            phaseSpaceState[count++] = p.v.y;
        }
    }
    
    /**
     * Gets the dimension of the phase space state
     * (particles * 2 dimensions * 2 for velocity and position)
     * @return dimension
     */
    public int getPhaseSpaceDim() {        

        return info.particles.size() * 4;
    }
    
    /**
     * Sets the phase space state of the particle system
     * @param phaseSpaceState
     */
    public void setPhaseSpace( double[] phaseSpaceState ) {
        int count = 0;
        for ( Particle particle : info.particles ) {
            if ( particle.pinned ) {
                count += 4;
            } else {
                particle.p.x = phaseSpaceState[count++];
                particle.p.y = phaseSpaceState[count++];
                particle.v.x = phaseSpaceState[count++];
                particle.v.y = phaseSpaceState[count++];
            }
        }
    }
    
    /**
     * Fixes positions and velocities after a step to deal with collisions 
     */
    public void postStepFix() {
        for ( Particle p : info.particles ) {
            if ( p.pinned ) {
                p.v.set(0,0);
            }
        }
        // do wall collisions
        double r = restitution.getValue();
        for ( Particle p : info.particles ) {
            if ( p.p.x <= 0 ) {
                p.p.x = 0;
                if ( p.v.x < 0 ) p.v.x = - p.v.x * r;
                if ( p.f.x < 0 ) p.f.x = 0;                
            }
            if ( p.p.x >= width ) {
                p.p.x = width;
                if (p.v.x > 0 ) p.v.x = - p.v.x * r;
                if (p.f.x > 0 ) p.f.x = 0;
            } 
            
            if ( p.p.y >= height ) {
                p.p.y = height;
                if ( p.v.y > 0 ) p.v.y = - p.v.y * r;
                if ( p.f.y > 0 ) p.f.y = 0;
            } 
            if ( p.p.y <= 0 ) {
                p.p.y = 0;
                if ( p.v.y < 0 ) p.v.y = - p.v.y * r;
                if ( p.f.y < 0 ) p.f.y = 0;
            }
        }
    }
    
    /** Elapsed simulation time */
    public double time = 0;

    /** The explicit integrator to use, if not performing backward Euler implicit integration */
    public Integrator integrator;
    
    public double[] state = new double[1];
    public double[] stateOut = new double[1];

    // these get created in init() and are probably useful for Backward Euler computations
    private ConjugateGradientMTJ CG;
    private DenseMatrix A;
    private DenseMatrix dfdx;
    private DenseMatrix dfdv;
    private DenseVector deltaxdot;
    private DenseVector b;
    private DenseVector f;
    private DenseVector xdot;
    
    /**
     * Initializes the system 
     * Allocates the arrays and vectors necessary for the solve of the full system
     */
    public void init() {
        int N = info.particles.size();
        // create matrix and vectors for solve
        CG = new ConjugateGradientMTJ(2*N);
        CG.setFilter(this);
        A = new DenseMatrix(2*N, 2*N);
        dfdx = new DenseMatrix(2*N, 2*N);
        dfdv = new DenseMatrix(2*N, 2*N);
        deltaxdot = new DenseVector(2*N);
        b = new DenseVector(2*N);
        f = new DenseVector(2*N);
        xdot = new DenseVector(2*N);
    }
    
    /**
     * Fills in the provided vector with the particle velocities.
     * @param xd
     */
    private void getVelocities(DenseVector xd) {
        for ( Particle p : info.particles ) {
            int j = p.index * 2;
            if( p.pinned ) {
                xd.set( j, 0 );
                xd.set( j+1, 0 );
            } else {
                xd.set( j, p.v.x );
                xd.set( j+1, p.v.y );
            }
        }       
    }

    /**
     * Sets the velocities of the particles given a vector
     * @param xd
     */
    private void setVelocities(DenseVector xd) {
        for ( Particle p : info.particles ) {
            int j = p.index * 2;
            if( p.pinned ) {
                p.v.set(0,0);
            } else {
                p.v.x = xd.get(j);
                p.v.y = xd.get(j+1);
            }
        }
    }
    
    /**
     *  Evaluates derivatives for ODE integration.
     * @param t time 
     * @param p phase space state
     * @param dpdt to be filled with the derivative
     */
    @Override
    public void derivs(double t, double[] p, double[] dpdt) {
        // set particle positions to given values
        setPhaseSpace( p );
        
        // TODO: Objective 2, for explicit integrators, compute forces, and accelerations, and set dpdt

        for(Particle particle: info.particles)
        {
            particle.clearForce();

            if(useGravity.getValue())
            {
                particle.f.y = gravity.getValue() * particle.mass;
            }

            particle.f.x -= viscousDamping.getValue() * particle.v.x;
            particle.f.y -= viscousDamping.getValue() * particle.v.y;
        }

        for(Spring s: info.springs)
        {
            // Applies the spring force by adding a force to each particle
            s.apply();
        }

        // overwrite dpdt
        int i = 0;
        for(Particle particle: info.particles)
        {
            dpdt[i++] = particle.v.x;
            dpdt[i++] = particle.v.y;
            // compute the accelerations from the forces
            dpdt[i++] = particle.f.x / particle.mass;
            dpdt[i++] = particle.f.y / particle.mass;
        }

    }
    
    /** Time in seconds that was necessary to advance the system */
    public double computeTime;
    
    /**
     * Advances the state of the system
     * @param elapsed
     */
    public void advanceTime( double elapsed ) {
        Spring.k = springStiffness.getValue();
        Spring.c = springDamping.getValue();
            
        int n = getPhaseSpaceDim();
        
        long now = System.nanoTime();        
        
        if ( explicit.getValue() ) {
            if ( n != state.length ) {
                state = new double[n];
                stateOut = new double[n];
            }
            // TODO: See explicit stepping here
            getPhaseSpace(state);         
            integrator.step( state, n, time, elapsed, stateOut, this);                
            setPhaseSpace(stateOut);
        } else {        
            if ( f == null || f.size() != n ) {
                init();
            }
            
            // TODO: Objective 8, your backward Euler implementation will go here!
            // Note that the init() method called above creates a bunch of very 
            // useful MTJ working variables for you, and the ConjugateGradientMTJ object.
            // Go look at that code now!
            
            
            
            
        }
        time = time + elapsed;
        postStepFix();
        computeTime = (System.nanoTime() - now) / 1e9;
    }
    
    @Override
    public void filter(Vector v) {
        for ( Particle p : info.particles ) {
            if ( !p.pinned ) continue;
            v.set( p.index*2+0, 0 );
            v.set( p.index*2+1, 0 );
        }
    }

    /**
     * Creates a new particle and adds it to the system
     * @param x
     * @param y
     * @param vx
     * @param vy
     * @return the new particle
     */
    public Particle createParticle( double x, double y, double vx, double vy ) {
        Particle p = new Particle( x, y, vx, vy );
        p.index = info.particles.size();
        info.particles.add( p );
        return p;
    }
    
    public void remove( Particle p ) {
    	for ( Spring s : p.springs ) {
    		Particle other = s.p1 == p ? s.p2 : s.p1; 
    		other.springs.remove( s );
            info.springs.remove( s );
    	}
    	p.springs.clear(); // not really necessary
        info.particles.remove( p );
    	// reset indices of each particle :(
    	for ( int i = 0 ; i < info.particles.size(); i++ ) {
            info.particles.get(i).index = i;
    	}
    }
    
    /**
     * Creates a new spring between two particles and adds it to the system.
     * @param p1
     * @param p2
     * @return the new spring
     */
    public Spring createSpring( Particle p1, Particle p2 ) {
        Spring s = new Spring( p1, p2 );
        info.springs.add( s );
        return s;
    }
    
    /**
     * Removes a spring between p1 and p2 if it exists, does nothing otherwise
     * @param p1
     * @param p2
     * @return true if the spring was found and removed
     */
    public boolean removeSpring( Particle p1, Particle p2 ) {
    	Spring found = null;
    	for ( Spring s : info.springs ) {
    		if ( ( s.p1 == p1 && s.p2 == p2 ) || ( s.p1 == p2 && s.p2 == p1 ) ) {
    			found = s;
    			break;
    		}
    	}
    	if ( found != null ) {
    		found.p1.springs.remove(found);
    		found.p2.springs.remove(found);
            info.springs.remove(found);
			return true;
    	}
    	return false;
    }
    
    @Override
    public void init(GLAutoDrawable drawable) {
        // do nothing
    }

    private int height;
    private int width;

    @Override
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();

        // update the width and the height for wall collisions
        height = drawable.getSurfaceHeight();
        width = drawable.getSurfaceWidth();
        
        gl.glPointSize( 10 );
        gl.glBegin( GL.GL_POINTS );
        for ( Particle p : info.particles ) {
            double alpha = 0.5;
            if ( p.pinned ) {
                gl.glColor4d( 1, 0, 0, alpha );
            } else {
                gl.glColor4d( p.color.x, p.color.y, p.color.z, alpha );
            }
            gl.glVertex2d( p.p.x, p.p.y );
        }
        gl.glEnd();
        
        gl.glColor4d(0,.5,.5,.5);
        gl.glLineWidth(2f);
        gl.glBegin( GL.GL_LINES );
        for (Spring s : info.springs) {
            gl.glVertex2d( s.p1.p.x, s.p1.p.y );
            gl.glVertex2d( s.p2.p.x, s.p2.p.y );
        }
        gl.glEnd();
    }
    
    public BooleanParameter useGravity = new BooleanParameter( "use gravity", true );
    public DoubleParameter gravity = new DoubleParameter( "gravity", 9.8, 0.01, 1000 );
    public DoubleParameter springStiffness = new DoubleParameter( "spring stiffness", 1, 0, 10000 );
    public DoubleParameter springDamping = new DoubleParameter( "spring damping", 0.1, 0, 50 );
    public DoubleParameter viscousDamping = new DoubleParameter( "viscous damping", 0.1, 0, 10 );
    public DoubleParameter restitution = new DoubleParameter( "r", 0, 0, 1 );
    public JTextArea comments = new JTextArea("enter comments in control panel");
    public IntParameter iterations = new IntParameter( "iterations", 100, 1, 100 );
    /** controls weather explicit or implicit integration is used */
    public BooleanParameter explicit = new BooleanParameter( "explicit", true );
    
    @Override
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.add( comments );
        vfp.add( useGravity.getControls() );
        vfp.add( gravity.getSliderControls(true) );
        vfp.add( springStiffness.getSliderControls(false) );
        vfp.add( springDamping.getSliderControls(false) );
        vfp.add( viscousDamping.getSliderControls(false) );
        vfp.add( restitution.getSliderControls(false) );
        vfp.add( iterations.getSliderControls() );
        vfp.add( explicit.getControls() );
        return vfp.getPanel();        
    }
    
    @Override
    public String toString() {
        String ret = "Marie Cornellier\n" +
                     comments.getText() + "\n" +
                     "particles = " + info.particles.size() + "\n";
        if ( explicit.getValue() ) {
            ret += "integrator = " + integrator.getName() + "\n";
        } else {
            ret += "integrator = Backward Euler\n";
        }
        ret += "k = " + springStiffness.getValue() + "\n" +
               "c = " + springDamping.getValue() + "\n" +
               "b = " + viscousDamping.getValue() +"\n" + 
               "time = " + time;
        return ret;
    }
    
}
