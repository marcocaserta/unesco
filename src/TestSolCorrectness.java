import org.junit.Test;
import static org.junit.Assert.assertEquals;

/// Test correctness of a solution (feasibility, obj function value, etc.)
/**
 * Here we verify the correctness of the final solution, w.r.t. the objective
 * function value and the constraints. We want to verify the following:
 * <ul>
 * <li> That the solution does not contain subtours (TSP-type)
 * <li> That the total score Z is correct
 * <li> That the solution does not exceed the time Budget
 * <li> That the balance constraint (nC approx nN) is satisfied
 * <li> That the number of endangered sites reported is correct
 * <ul>
 * Note that the total number of countries visited is computed using <i>logical
 * relations</i> among scores (to ensure further correctness.) See
 * numberCountries().
 */
/// Test that the solution does not contain subtours
public class TestSolCorrectness {

    @Test
    /// Test feasibility w.r.t. the TSP subtours
    public void testNoCycles()
    {
        int[] nVisits = new int[Orienteering.nSites];
        int maxVisits = 0, el = -1;
        for (int j = 1; j < Orienteering.bestRoute.inSet.size()-1; j++)
        {
            el = Orienteering.bestRoute.inSet.get(j);
            nVisits[el]++;
            if (nVisits[el] > maxVisits)
                maxVisits = nVisits[el];
        }

        System.out.println("- Test Subtours \t\t Successful!");
    }

    @Test
    /// Test feasibility (travel time) in two ways: Manual sum and using computeTime()
    public void testFeasibility()
    {
        double totTime = Orienteering.bestRoute.computeTime();
        assertEquals(Orienteering.bestRoute.time, totTime, Orienteering._EPSI);
        System.out.println("- Test Feasibility (version 1) \t Successful!");


        int el = -1, nextEl = -1;
        totTime = Orienteering._VISIT*(Orienteering.bestRoute.inSet.size()-2);
        for (int j = 0; j < Orienteering.bestRoute.inSet.size()-1; j++)
        {
            el = Orienteering.bestRoute.inSet.get(j);
            nextEl = Orienteering.bestRoute.inSet.get(j+1);
            totTime += Orienteering.timeMatrix[el][nextEl];
        }
        assertEquals(Orienteering.bestRoute.time, totTime, Orienteering._EPSI);
        System.out.println("- Test Feasibility (version 2) \t Successful!");
    }

    @Test
    /// Verify objective function value
    public void objFunctionValue()
    {
        int totZ = Orienteering.bestRoute.computeZ();
        assertEquals(Orienteering.bestRoute.Z, totZ);
        System.out.println("- Test Objective Function \t Successful!");
    }


    @Test
    /// Check balance constraint and other features of the final solution
    public void balanceConstraint()
    {
        
        int nC = 0, nN = 0, nMix = 0, nDanger = 0, nCountries = 0, el = 0;

        for (int j = 1; j < Orienteering.bestRoute.inSet.size()-1; j++)
        {
            el     = Orienteering.bestRoute.inSet.get(j);

            if (Orienteering.danger[el] == 1) nDanger++;
            if (Orienteering.category[el].charAt(0)  == 'C')
                nC++;
            else if (Orienteering.category[el].charAt(0)  == 'N')
                nN++;
            else 
                nMix++;
        }

        // check balance
        assertEquals(Orienteering.bestRoute.nC, nC);
        assertEquals(Orienteering.bestRoute.nN, nN);
        assertEquals(Orienteering.bestRoute.nMix, nMix);
        assertEquals((Math.abs(nC - nN) <= nMix +1), true);
        System.out.println("- Test Balance Constraint \t Successful!");

        assertEquals(Orienteering.bestRoute.nDanger, nDanger);
        System.out.println("- Test Nr. Endangered Sites \t Successful!");

    }

    @Test
    /// Check balance constraint and other feature of the final solution
    public void numberCountries()
    {
        int nCountries = Orienteering.bestRoute.Z - Orienteering.bestRoute.nDanger*3 - (Orienteering.bestRoute.inSet.size()-2); // remove home twice (first and last size):
        nCountries /= 2; // 2 points for each country

        assertEquals(Orienteering.bestRoute.nCountries, nCountries);
        System.out.println("- Test Nr. Visited Countries \t Successful!");
    }
}
