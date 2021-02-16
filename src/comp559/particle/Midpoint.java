package comp559.particle;

public class Midpoint implements Integrator {

    @Override
    public String getName() {
        return "midpoint";
    }

    private double[] midpoints;
    private double[] dpdt;
    
    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
        // TODO: Objective 4, implement midpoint method

    	// You will probably want a temporary array in this method and perhaps
    	// multiple temporary arrays in other higher order explicit integrators.
    	// Avoid thrashing memory by reallocating only when necessary.
    	if ( midpoints == null || midpoints.length != n ) {
            midpoints = new double[n];
    	}
        if ( dpdt == null || dpdt.length != n ) {
            dpdt = new double[n];
        }

        derivs.derivs(t, p, dpdt); // get the slope at current phase state and current time

        for(int i = 0; i < n; i++)
        {
            midpoints[i] = p[i] + (h/2.0)*dpdt[i]; // step halfway
        }

        derivs.derivs(t+h/2.0, midpoints, dpdt); // get new slope at half time step later, with the midpoint phase state

        for(int i = 0; i < n; i++)
        {
            pout[i] = p[i] + h* dpdt[i];
        }
    	
    }

}
