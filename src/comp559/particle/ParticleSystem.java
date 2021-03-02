package comp559.particle;

import java.util.Iterator;
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
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.SceneGraphNode;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;

import java.io.*;

import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

/**
 * Implementation of a simple particle system
 * @author kry
 */
public class ParticleSystem implements SceneGraphNode, Function, Filter {

    public List<Particle> particles;
    public List<Spring> springs;
    private Matrix K;
    private DenseVector x0;
    private DenseVector v0;
    private DenseVector xNew;
    private DenseVector vNew;
    private DenseVector tmp;
    
    /**
     * Creates an empty particle system
     */
    public ParticleSystem() {
        particles = new LinkedList<Particle>();
        springs = new LinkedList<Spring>();
    }

    /**
     * Create a saved test system from the file
     * @param file_location
     */
    public void loadSystem( String file_location ) {
        try {
            File excel = new File(file_location);
            FileInputStream fis = new FileInputStream(excel);
            XSSFWorkbook book = new XSSFWorkbook(fis);
            XSSFSheet sheet = book.getSheetAt(0);

            Iterator<Row> itr = sheet.iterator();
            Row row = itr.next();
            Iterator<Cell> cellIterator = row.cellIterator();
            Cell cell = cellIterator.next();
            cell = cellIterator.next();
            int num_particles = (int) cell.getNumericCellValue();
            row = itr.next();

            Particle[] loaded_particles = new Particle[num_particles];

            // Iterating over Excel file in Java
            int i = 0;
            while (itr.hasNext()) {
                row = itr.next();

                // Iterating over each column of Excel file
                cellIterator = row.cellIterator();
                cell = cellIterator.next();
                cell = cellIterator.next();
                double pos_x = cell.getNumericCellValue();
                cell = cellIterator.next();
                double pos_y = cell.getNumericCellValue();
                cell = cellIterator.next();
                double mass = cell.getNumericCellValue();

                Particle p = new Particle(pos_x, pos_y, 0, 0, mass);
                p.index = i;
                particles.add( p );

                loaded_particles[i++] = p;
            }

            XSSFSheet sheet2 = book.getSheetAt(1);

            Iterator<Row> itr2 = sheet2.iterator();
            row = itr2.next();
            cellIterator = row.cellIterator();
            cell = cellIterator.next();
            cell = cellIterator.next();
            int num_springs = (int) cell.getNumericCellValue();
            row = itr2.next();

            while (itr2.hasNext()) {
                row = itr2.next();

                // Iterating over each column of Excel file
                cellIterator = row.cellIterator();
                cell = cellIterator.next();
                cell = cellIterator.next();
                int p1_index = (int) cell.getNumericCellValue();
                cell = cellIterator.next();
                int p2_index = (int) cell.getNumericCellValue();
                cell = cellIterator.next();
                double k_value = cell.getNumericCellValue();

                Spring s = new Spring(loaded_particles[p1_index], loaded_particles[p2_index]);
                springs.add(s);
            }
        } catch (FileNotFoundException fe) {
            fe.printStackTrace();
        } catch (IOException ie) {
            ie.printStackTrace();
        }

    }

    /**
     * Saves current particle system to a file that can be loaded later
     * @param filename
     */
    public void saveSystem( String filename ) {
        XSSFWorkbook workbook = new XSSFWorkbook();
        XSSFSheet sheet = workbook.createSheet("Particles");
        XSSFSheet sheet2 = workbook.createSheet("Springs");

        int rowCount = 0;
        Row row = sheet.createRow(0);
        Cell cell = row.createCell(0);
        cell.setCellValue("Particles");
        cell = row.createCell(1);
        cell.setCellValue(particles.size());

        row = sheet.createRow(++rowCount);
        cell = row.createCell(0);
        cell.setCellValue("Index");
        cell = row.createCell(1);
        cell.setCellValue("x");
        cell = row.createCell(2);
        cell.setCellValue("y");
        cell = row.createCell(3);
        cell.setCellValue("mass");

        for (Particle p : particles) {
            row = sheet.createRow(++rowCount);

            cell = row.createCell(0);
            cell.setCellValue(p.index);
            cell = row.createCell(1);
            cell.setCellValue(p.p0.x);
            cell = row.createCell(2);
            cell.setCellValue(p.p0.y);
            cell = row.createCell(3);
            cell.setCellValue(p.mass);

        }

        rowCount = 0;
        row = sheet2.createRow(0);
        cell = row.createCell(0);
        cell.setCellValue("Springs");
        cell = row.createCell(1);
        cell.setCellValue(springs.size());

        row = sheet2.createRow(++rowCount);
        cell = row.createCell(0);
        cell.setCellValue("Index");
        cell = row.createCell(1);
        cell.setCellValue("p1");
        cell = row.createCell(2);
        cell.setCellValue("p2");
        cell = row.createCell(3);
        cell.setCellValue("k");

        for (Spring s : springs) {
            row = sheet2.createRow(++rowCount);

            cell = row.createCell(0);
            cell.setCellValue(rowCount - 2);
            cell = row.createCell(1);
            cell.setCellValue(s.p1.index);
            cell = row.createCell(2);
            cell.setCellValue(s.p2.index);
            cell = row.createCell(3);
            cell.setCellValue(s.k);

        }

        try {
            File file = new File("./savedSystems/" + filename);
            FileOutputStream outputStream = new FileOutputStream(file);
            workbook.write(outputStream);
        } catch (FileNotFoundException fe) {
            fe.printStackTrace();
        } catch (IOException ie) {
            ie.printStackTrace();
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
            p1.index = particles.size();
            particles.add( p1 );
            p2 = new Particle( p.x + d.y, p.y - d.x, 0, 0 );
            p2.index = particles.size();
            particles.add( p2 );
            springs.add( new Spring ( p1, p2 ) );
            p1.pinned = true;
            p2.pinned = true;            
            p.add( d );
            p.add( d );                    
            int N = 10;
            for (int i = 1; i < N; i++ ) {                
                //d.set( 20*Math.cos(i*Math.PI/N), 20*Math.sin(i*Math.PI/N) );                
                p3 = new Particle( p.x - d.y, p.y + d.x, 0, 0 );
                p3.index = particles.size();
                particles.add( p3 );
                p4 = new Particle( p.x + d.y, p.y - d.x, 0, 0 );
                p4.index = particles.size();
                particles.add( p4 );
                springs.add( new Spring ( p3, p1 ) );
                springs.add( new Spring ( p3, p2 ) );
                springs.add( new Spring ( p4, p1 ) );
                springs.add( new Spring ( p4, p2 ) );
                springs.add( new Spring ( p4, p3 ) );
                p1 = p3;
                p2 = p4;                
                p.add( d );
                p.add( d );            
            }
        } else if ( which == 2) {
            Particle p1 = new Particle( 320, 100, 0, 0 );
            p1.index = particles.size();
            particles.add( p1 );
            Particle p2 = new Particle( 320, 200, 0, 0 );
            p2.index = particles.size();
            particles.add( p2 );
            p1.pinned = true;
            springs.add( new Spring( p1, p2 ) );
        } else if ( which == 3 ) {
            int ypos = 100;
            Particle p0 = null;
            Particle p1, p2;
            p1 = new Particle( 320, ypos, 0, 0 );
            p1.index = particles.size();
            p1.pinned = true;
            particles.add( p1 );
            int N = 10;
            for ( int i = 0; i < N; i++ ) {
                ypos += 20;
                p2 = new Particle( 320, ypos, 0, 0 );
                p2.index = particles.size();
                particles.add( p2 );
                springs.add( new Spring( p1, p2 ) );
                // Hum.. this is not great in comparison to a proper bending energy...
                // use Maple to generate some code though, as it is painful to write by hand! :(
                if ( p0 != null ) springs.add( new Spring( p2, p0 ) );
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
        return particles;
    }
    
    /**
     * Gets the springs in the system
     * @return the spring list
     */
    public List<Spring> getSprings() {
    	return springs;
    }
    
    /**
     * Resets the positions of all particles
     */
    public void resetParticles() {
        for ( Particle p : particles ) {
            p.reset();
        }
        time = 0;
    }
    
    /**
     * Deletes all particles
     */
    public void clearParticles() {
        particles.clear();
        springs.clear();
    }
    
    /**
     * Gets the phase space state of the particle system
     * @param phaseSpaceState
     */
    public void getPhaseSpace( double[] phaseSpaceState ) {
        int count = 0;
        for ( Particle p : particles ) {
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

        return particles.size() * 4;
    }
    
    /**
     * Sets the phase space state of the particle system
     * @param phaseSpaceState
     */
    public void setPhaseSpace( double[] phaseSpaceState ) {
        int count = 0;
        for ( Particle particle : particles ) {
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
        for ( Particle p : particles ) {
            if ( p.pinned ) {
                p.v.set(0,0);
            }
        }
        // do wall collisions
        double r = restitution.getValue();
        // TODO: dpi stuff
        for ( Particle p : particles ) {
            if ( p.p.x - (5 * p.mass / 2) <= 0 ) {
                p.p.x = (5 * p.mass / 2);
                if ( p.v.x < 0 ) p.v.x = - p.v.x * r;
                if ( p.f.x < 0 ) p.f.x = 0;                
            }
            if ( p.p.x >= width - (5 * p.mass / 2) ) {
                p.p.x = width - (5 * p.mass / 2);
                if (p.v.x > 0 ) p.v.x = - p.v.x * r;
                if (p.f.x > 0 ) p.f.x = 0;
            } 
            
            if ( p.p.y >= height - (5 * p.mass / 2) ) {
                p.p.y = height - (5 * p.mass / 2);
                if ( p.v.y > 0 ) p.v.y = - p.v.y * r;
                if ( p.f.y > 0 ) p.f.y = 0;
            } 
            if ( p.p.y - (5 * p.mass / 2) <= 0 ) {
                p.p.y = (5 * p.mass / 2);
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
    private FlexCompRowMatrix I;
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
        int N = particles.size();
        // create matrix and vectors for solve
        CG = new ConjugateGradientMTJ(2*N);
        CG.setFilter(this);
        A = new DenseMatrix(2*N, 2*N);
        I = new FlexCompRowMatrix(2*N, 2*N);
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
        for ( Particle p : particles ) {
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
        for ( Particle p : particles ) {
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

        for(Particle particle: particles)
        {
            particle.clearForce();

            // apply force due to gravity on each particle
            if(useGravity.getValue())
            {
                particle.f.y = gravity.getValue() * particle.mass;
            }

            // viscous damping is a force directly applied by the environment
            particle.f.x -= viscousDamping.getValue() * particle.v.x;
            particle.f.y -= viscousDamping.getValue() * particle.v.y;
        }

        for(Spring s: springs)
        {
            // Applies the spring force by adding an equal force to each particle
            s.apply();
        }

        // overwrite dpdt
        int i = 0;
        for(Particle particle: particles)
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
        Spring.default_k = springStiffness.getValue();
        Spring.default_c = springDamping.getValue();
            
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
            if ( n != state.length ) {
                state = new double[n];
                stateOut = new double[n];
            }
            
            // TODO: Objective 8, your backward Euler implementation will go here!
            // Note that the init() method called above creates a bunch of very 
            // useful MTJ working variables for you, and the ConjugateGradientMTJ object.
            // Go look at that code now!

            getPhaseSpace(state);
            backwardEuler(state, n, elapsed, stateOut);
            setPhaseSpace(stateOut);
            
        }
        time = time + elapsed;
        postStepFix();
        computeTime = (System.nanoTime() - now) / 1e9;
    }
    
    @Override
    public void filter(Vector v) {
        for ( Particle p : particles ) {
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
        Particle p = new Particle( x, y, vx, vy, mass.getValue());
        p.index = particles.size();
        particles.add( p );
        return p;
    }

    /**
     * Creates an imaginary mouse particle and adds it to the system
     * @param x
     * @param y
     * @param vx
     * @param vy
     * @return the new particle
     */
    public Particle createMouseParticle( double x, double y, double vx, double vy ) {
        Particle p = new Particle( x, y, vx, vy,1,false );
        p.index = particles.size();
        particles.add( p );
        return p;
    }
    
    public void remove( Particle p ) {
    	for ( Spring s : p.springs ) {
    		Particle other = s.p1 == p ? s.p2 : s.p1; 
    		other.springs.remove( s );
            springs.remove( s );
    	}
    	p.springs.clear(); // not really necessary
        particles.remove( p );
    	// reset indices of each particle :(
    	for ( int i = 0 ; i < particles.size(); i++ ) {
            particles.get(i).index = i;
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
        springs.add( s );
        return s;
    }

    /**
     * Creates the mouse spring between two particles and adds it to the system.
     * @param p1
     * @param p2
     * @return the new spring
     */
    public Spring createMouseSpring( Particle p1, Particle p2, double k, double c, double l0) {
        Spring s = new Spring( p1, p2, k, c, l0, false);
        springs.add( s );
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
    	for ( Spring s : springs ) {
    		if ( ( s.p1 == p1 && s.p2 == p2 ) || ( s.p1 == p2 && s.p2 == p1 ) ) {
    			found = s;
    			break;
    		}
    	}
    	if ( found != null ) {
    		found.p1.springs.remove(found);
    		found.p2.springs.remove(found);
            springs.remove(found);
			return true;
    	}
    	return false;
    }


    private void backwardEuler(double[] state, int n, double h, double[] stateOut)
    {
        int N = particles.size();
        x0 = new DenseVector(2*N);
        v0 = new DenseVector(2*N);
        xNew = new DenseVector(2*N);
        vNew = new DenseVector(2*N);
        tmp = new DenseVector(2*N);

        for(int i = 0; i < 2*N; i++)
        {
            I.set(i, i, 1.0);
        }

        for(int i = 0; i < state.length; i+=4)
        {
            int j = i/2;
            x0.set(j, state[i]);
            x0.set(j+1, state[i+1]);
            v0.set(j, state[i+2]);
            v0.set(j+1, state[i+3]);
        }

        for(Spring s: springs)
        {
            // apply all the spring forces
            s.addForce(f);
            s.addDfdx(dfdx);
            s.addDfdv(dfdv);
        }

        for(Particle p: particles)
        {
            f.add(p.index*2, -viscousDamping.getValue() * p.v.x);
            f.add(p.index*2+1, -viscousDamping.getValue() * p.v.y);
            f.add(p.index*2+1, gravity.getValue() * p.mass);
        }

        dfdx.mult(v0, tmp);
        A.set(I.add(dfdv.scale(-h)).add(dfdx.scale(-h*h)));
        b.set(f.add(tmp.scale(h)).scale(h));

        CG.solve(A, b, deltaxdot, iterations.getValue());

        vNew.set(v0.add(deltaxdot));
        tmp.set(vNew);
        xNew.set(x0.add(tmp.scale(h)));

        for(int i = 0; i < stateOut.length; i+=4)
        {
            int j = i/2;
            stateOut[i] = xNew.get(j);
            stateOut[i+1] = xNew.get(j+1);
            stateOut[i+2] = vNew.get(j);
            stateOut[i+3] = vNew.get(j+1);
        }
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
        
//        gl.glPointSize( 10 );

        for ( Particle p : particles ) {
            double alpha = 0.5;
            if (p.visible) {
                Vector2d v = new Vector2d(p.v.x, p.v.y);
                v.absolute();
                if ( p.pinned ) {
                    gl.glColor4d( 1, 1, 1, alpha );
//                    gl.glColor4d( v.x, v.y, 1, alpha );

                } else {
//                    gl.glColor4d( p.color.x, p.color.y, p.color.z, alpha );
                    gl.glColor4d( v.x, v.y, 1, alpha );
                }
                float sz = (float) (p.mass * 10);
                gl.glPointSize( sz );
                gl.glBegin( GL.GL_POINTS );
                gl.glVertex2d( p.p.x, p.p.y );
                gl.glEnd();
            }
        }

        gl.glColor4d(0,.5,.5,.5);
        gl.glLineWidth(2f);
        gl.glBegin( GL.GL_LINES );
        for (Spring s : springs) {
            if (s.visible) {
                gl.glVertex2d( s.p1.p.x, s.p1.p.y );
                gl.glVertex2d( s.p2.p.x, s.p2.p.y );
            }
        }
        gl.glEnd();
    }

    public DoubleParameter mass = new DoubleParameter( "mass", 1, 0.1, 100 );
    public BooleanParameter useGravity = new BooleanParameter( "use gravity", true );
    public DoubleParameter gravity = new DoubleParameter( "gravity", 9.8, 0.01, 1000 );
    public DoubleParameter springStiffness = new DoubleParameter( "spring stiffness", 1, 0, 10000 );
    public DoubleParameter springDamping = new DoubleParameter( "spring damping", 0.1, 0, 50 );
    public DoubleParameter viscousDamping = new DoubleParameter( "viscous damping", 0.1, 0, 10 );
    public DoubleParameter restitution = new DoubleParameter( "r", 0.5, 0, 1 );
    public JTextArea comments = new JTextArea("enter comments in control panel");
    public IntParameter iterations = new IntParameter( "iterations", 100, 1, 100 );
    /** controls weather explicit or implicit integration is used */
    public BooleanParameter explicit = new BooleanParameter( "explicit", true );
    
    @Override
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.add( comments );
        vfp.add( useGravity.getControls() );
        vfp.add( mass.getSliderControls(true) );
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
                     "particles = " + particles.size() + "\n";
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
