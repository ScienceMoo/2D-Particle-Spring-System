package comp559.particle;

public class ModifiedMidpoint implements Integrator {

    @Override
    public String getName() {
        return "modified midpoint";
    }

    private double[] midpoints;
    private double[] dpdt;

    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
    	// TODO: Objective 5, implemement the modified midpoint (2/3) method.
    	// see also efficient memory management suggestion in provided code for the Midpoint method.
        if ( midpoints == null || midpoints.length != n ) {
            midpoints = new double[n];
        }
        if ( dpdt == null || dpdt.length != n ) {
            dpdt = new double[n];
        }

        derivs.derivs(t, p, dpdt);

        for(int i = 0; i < n; i++)
        {
            midpoints[i] = p[i] + h*dpdt[i]*(2.0/3.0);
        }

        derivs.derivs(t+h*(2.0/3.0), midpoints, dpdt);

        for(int i = 0; i < n; i++)
        {
            pout[i] = p[i] + h*dpdt[i];
        }
    }

}
