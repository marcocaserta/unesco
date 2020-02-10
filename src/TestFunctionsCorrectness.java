import org.junit.Test;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import java.util.*;

/// Test correctness of a solution (feasibility, obj function value, etc.)
/**
 * Here we verify the correctness of the final solution, w.r.t. the objective
 * function value and the constraints. We want to verify the following:
 * <ul>
 * <li> That the total score Z is correct
 * <li> That the solution does not exceed the time Budget
 * <li> That the balance constraint (nC approx nN) is satisfied
 * <li> That the number of endangered sites reported is correct
 * <ul>
 * Note that the total number of countries visited is computed using logical
 * relations among scores (to ensure further correctness.) See
 * numberCountries().
 */
public class TestFunctionsCorrectness {

    @Test
    /// Test feasibility in two ways: Manual sum and using computeTime()
    /** We test the computation of a simple solution. Consider the following:
     * \f[ \textbf{x} = (1, 2, 3, 1)\f]
     * where 1 is home here (we do not use 0, since it can be changed by the
     * user.) The longitude and latitute of each site are given below:
     * <ul>
     * <li> Site 1 :: 67.82525,	34.84694, country = af, danger = 1  
     * <li> Site 2 :: 64.5158888889,	34.3964166667, country=af, danger = 1
     * <li> Site 3 :: 20.13333333,	40.06944444, country=al, danger = 0
     * </ul>
     * The objective function is computed as follows:
     * <ul>
     * <li> nr. sites = 2 (sites 2 and 3, since site 1 is home)
     * <li> nr. countries = 2 (af, al). Country of home is not counted
     * <li> nr. endangered = 1
     * </ul>
     * Therefore, we get:
     * \f[ z = 2 + 2\times 2 +3\times 1 = 9 \f]
     * In terms of distance, we computed the overall tour length using 
     * @see <a href="https://andrew.hedges.name/experiments/haversine/">https://andrew.hedges.name/experiments/haversine/</a>
     * and we obtained a total distance of 8477.141 km. We transform this
     * distance in time (using 80km/h as traveling speed) and we add 2*_VISIT
     * time for visiting (6 hours in each site.) We then compare the manual
     * computaton with the one provided by the function computeTime();
     *
     * The manual computation of time is as follows:
     * \f[ t(\textbf{x)= \frac{8477.141}{s} + 2\times 360 = 7077.85575
     * \f]
     * where \f$s\f$ is the speed in km/min. equal to \f$s = 80/60\f$, and 360
     * indicates the time spent visiting each site (6 hours per site.)
     *
     * The computation obtained using timeMatrix is described in the function
     * computeTime() and sums up the different legs plus the visiting time.
     *
     * The test is aimed at ensuring that the maximum gap between the manual
     * computation and the value provided by computeTime() is no more than 1%.
     *
     * */
    public void testFeasibilityFunction()
    {
        double maxGap = 0.01;

        ArrayList<Integer> sol = new ArrayList<Integer>();

        Solution testSol = new Solution(Orienteering.home, Orienteering.nSites);

        testSol.inSet.clear(); // remove 0 (added by default)
        testSol.inSet.add(1); // 67.82525,	34.84694  (1 is home here, not counted)
        testSol.inSet.add(2); // 64.5158888889,	34.3964166667 (af, danger = 1)
        testSol.inSet.add(3); // 20.13333333,	40.06944444  (albania, danger = 0)
        testSol.inSet.add(1); // 67.82525,	34.84694 (1 is home, skipped)

        int myZ = 2 + 2*2 + 3;
        assertEquals(myZ, testSol.computeZ());
        System.out.println("- Test Obj Function Computation \t Successful!");

        double totTime = 8477.141 / Orienteering._SPEED + 2*Orienteering._VISIT;
        double gap = (totTime - testSol.computeTime())/totTime;

        assertEquals(0.0, gap, maxGap);
        System.out.println("- Test Budget Function Computation \t Successful!");


    }


}
