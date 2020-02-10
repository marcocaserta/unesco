/**
 * \mainpage Orienteering : The Unesco Challenge
 * \authors Marco Caserta
 *
 * This project can be compiled using:
 * ~~~
 * javac -d . -sourcepath src src/Orienteering.java
 * ~~~
 * and executed using (with default options - use flag -h for help):
 * ~~~
 * java Orienteering
 * ~~~
 */
import java.util.*;
import java.util.Scanner;

import org.junit.runner.JUnitCore;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

/// Project Main Class
/**
 * This code solves a variant of the Orienteering problem. A thorough
 * presentation of the model is found in the document `summary.pdf` under the
 * folder docs/summary. please, have a look at that document to find a
 * description of the algorithm.
 *
 *
 *
 * This class implements the four heuristics of the Iterative Local Search
 * algorithm, namely:
 * <ul>
 * <li> the insertionHeuristic()
 * <li> the twoOptHeuristic()
 * <li> the repairHeuristic()
 * <li> the randomRemoveOne()
 *
 */
public class Orienteering
{
    public static Random random; //!< random number generator

    static final int    _INFTY       = Integer.MAX_VALUE; //!< Infinity for integers
    static final int    _NMAX        = 10; //!< Nr. of cycles without improvements
    static final int    _TOTMAX      = 1000; //!< Maximum number of cycles overall
    static final double _MAX_PERTURB = 0.5; //!< Max % of perturbation
    static final double _EPSI        = 1e-06; //!< Epsilon value (almost zero)
    static final double _BUDGET      = 21*24*60;//!<  budget in minutes
    static final double _VISIT       = 6*60;  /* time spent in each site */
    static final double _SPEED       = 8.0/6.0; //!< Km per minute
    public static double latitude    = 51.533241; //!< Latitude Home point
    public static double longitude   = -0.104839;; //!< Longitude Home point
    public static int    stochastic  = 0; //!< set to 1 if stochastic, 0 if determ

    public static int home = 0; //!< default id of HOME node
    public static int nSites = 0; //<! total number of sites in the dataset


    // store the attributes of each site (populated from csv file in ReadCSV)
    public static ArrayList<Record> DataFrame = new ArrayList<Record>();

    public static double[] lat; //!< latitude of each site
    public static double[] lon; //!< longitude of each site
    public static int[] danger; //!< whether a site is endangered (1) or not (0)
    public static String[] country; //!< country iso code of each site
    public static String[] category; //!< category (C, N, or C/N) of each site

    public static double[][] timeMatrix; //!< nSites x nSites travel time matrix
    // note that these are travel times, not km

    // store best solution found
    public static Solution bestRoute = new Solution(home, nSites);


    /// Main Algorithm
    /** This algorithm implements an Iterative Local Search (ILS) heuristic scheme
     * for a variant of the Orienteering problem.
     */
    public static void main(String[] args) 
    {
        ParseCmdLine parser = new ParseCmdLine(); // command line parser
        try{ 
            parser.readCommandLine(args);
        }
        catch (IllegalArgumentException ex)
        {
            System.out.println(ex);
            parser.usage();
            return;
        }

        // set deterministic or stochastic version of the algorithm
        if (stochastic==0)
        {
            long seed = 27;
            random = new Random(seed);
        }
        else
        {
            random = new Random();
        }

        printParameters(); // print algorithmic parameters to screen

        ReadCSV reader = new ReadCSV(); // read data from csv file from disk
        reader.readData();

        // Initialization of data structures
        // a "cheaper" copy of the best solution (we only save list of sites)
        ArrayList<Integer> bestSet = new ArrayList<Integer>();
        int bestZ  = 0; // obj function of best solution
        nSites     = DataFrame.size()+1; // including HOME
        lon        = new double[nSites];
        lat        = new double[nSites];
        danger     = new int[DataFrame.size()+1];
        country    = new String[DataFrame.size()+1];
        category   = new String[DataFrame.size()+1];
        timeMatrix = new double[nSites][nSites]; // symmetric, even if triangular

        // save Home node as node 0
        lon[0]      = longitude;
        lat[0]      = latitude;
        danger[0]   = -1;
        country[0]  = " ";
        category[0] = " ";
        for (int i = 0; i < DataFrame.size(); i++)
        {
            lon[i+1]      = DataFrame.get(i).longitude;
            lat[i+1]      = DataFrame.get(i).latitude;
            danger[i+1]   = DataFrame.get(i).danger;
            country[i+1]  = DataFrame.get(i).iso_code;
            category[i+1] = DataFrame.get(i).category_short;
        }

        // compute DISTANCE MATRIX using HaversineDistance
        // NOTE: We convert distances in travel times (in minutes)
        computeDistances(timeMatrix);

        // UNITTEST # 1: run tests to verify correctness of functions
        Result result = JUnitCore.runClasses(TestFunctionsCorrectness.class);
        for (Failure failure : result.getFailures()) {
            System.out.println(failure.toString());
        }
        //System.out.println("");

        //----------------------------
        // Algorithm cycle starts here
        //----------------------------
        Solution route    = new Solution(home, nSites);
        boolean stopping  = false;
        int noImprovement = 0; //!< count nr iterations without improvements
        int nCycle        = 0; //!< iteration counter
        // stopping criteria: max cycles without improvement and overall max
        while ((noImprovement < _NMAX) && (nCycle < _TOTMAX))
        {
            nCycle++;
            insertionHeuristic(route); // greedy construction of a solution

            // restore feasibility w.r.t. number of C and N
            stopping = false;
            while (!stopping)
                stopping = repairHeuristic(route);

            // 2-Opt heuristic (ILS fashion)
            stopping = false;
            while (!stopping)
                stopping = twoOptHeuristic(route);

            route.computeSolScores(); // recomputing counting values
            // restore feasibility w.r.t. number of C and N
            stopping = false;
            while (!stopping)
                stopping = repairHeuristic(route);

            if (route.Z > bestZ) // update best solution found so far
            {
                bestZ         = route.Z;
                bestSet       = new ArrayList<Integer>(route.inSet);
                noImprovement = 0;
                System.out.printf("[Iter %3d] Best Solution :: Z = %5d (%3d sites)\n", nCycle, route.Z, route.inSet.size()-2);
            }
            else 
                noImprovement++; // if no improvement, increase counter

            // random perturbation of a % of the current solution
            routePerturbation(route, _MAX_PERTURB);
            route.computeSolScores();
        } // cycle until _NMAX cycles without improvements are reached

        // management of the best found solution
        bestRoute.inSet = new ArrayList<Integer>(bestSet);
        bestRoute.computeSolScores();
        printFinalRoute(bestRoute);


      // UNITTEST # 2: check correctness of the solution
      result = JUnitCore.runClasses(TestSolCorrectness.class);
      for (Failure failure : result.getFailures()) {
         System.out.println(failure.toString());
      }

      } 
    /************************ main program ******************************/
    /// END main program
    /************************ main program ******************************/

    /// Summary of algorithmic parameters.
    static void printParameters()
    {
        System.out.printf("\nUsing Algorithm Version : ");
        if (stochastic==0)
            System.out.println("** DETERMINISTIC **");
        else
            System.out.println("** STOCHASTIC **");
        System.out.println("");
        System.out.printf("Max Nr. Cycles \t \t = %8d\n", _NMAX);
        System.out.printf("Max Perturbation \t = %8.2f\n", _MAX_PERTURB);

        System.out.println("Starting Point ::");
        System.out.printf("\t Longitude \t = %8.4f\n", longitude);
        System.out.printf("\t Latitude  \t = %8.4f\n", latitude);
        System.out.printf("\n");
    }

    /// Perturb a percentage of the overall route.
    /**
     */
    static void routePerturbation(Solution route, double percent)
    {
        int nPerturb = (int)Math.ceil(percent*(double)route.inSet.size());
        int pos = -1;
        int nEls =  route.inSet.size();
        for (int t = 0; t < nPerturb; t++)
            {
                pos = 1 + random.nextInt(nEls-2);
                nEls--; // the solution vector shrinks by one at each iteration
                int el = route.inSet.get(pos);
                route.inSet.remove(pos);
                route.outSet.add(el);
            }
        route.computeSolScores();
    }

    /// Print best solution to screen.
    /**
     * Note: travelMatrix contains travel times, in minutes. We reconvert
     * these values back to km.
     * */
    static void printFinalRoute(Solution route)
    {
        System.out.println("\n\n----------------------");
        System.out.println(    "<< Best Route Found >> ");
        System.out.println(    "----------------------");
        System.out.println("Starting Point ::");
        System.out.printf("\t Longitude = %.4f\n", longitude);
        System.out.printf("\t Latitude  = %.4f\n", latitude);
        System.out.printf("\n");

        // REM.: Zero is Home. Therefore, when accessing dataFrame, we need
        // reduce the index by 1 (w.r.t. the value in 'route')
        int el, nEl, nElDF;
        String prefix, site, category, country;
        double totKm = 0.0;
        for (int j = 1; j < route.inSet.size()-1; j++)
        {
            el = route.inSet.get(j-1);
            nEl = route.inSet.get(j);
            nElDF = nEl - 1;
            totKm += timeMatrix[el][nEl]*_SPEED;
            site = DataFrame.get(nElDF).name_en;
            category = DataFrame.get(nElDF).category_short;
            country  = DataFrame.get(nElDF).states_name_en;
            if (DataFrame.get(nElDF).danger == 1) 
                prefix = "*";
            else
                prefix = " ";
            System.out.printf("%3d - [%5.0f km]  %s%s (%s) -- %s\n",j, _SPEED*timeMatrix[el][nEl], prefix, site, category, country);
        }

        System.out.printf("    - [%5.0f km]   Back Home.\n", _SPEED*timeMatrix[route.inSet.get(route.inSet.size()-2)][0]);
        totKm += _SPEED*timeMatrix[route.inSet.get(route.inSet.size()-2)][0];

        route.totKm = totKm;


        route.printSol();
        assert (Math.abs(route.Z - route.computeZ()) < _EPSI && 
                Math.abs(route.time - route.computeTime()) < _EPSI);

    }

    /// Compute distance matrix using haversine function.
    /** This function return the timeMatrix matrix. The matrix is triangular.
     * However, we keep here the full matrix, for ease of access. 
     *
     * NOTE: This is not a distance matrix. The values are travel times across
     * points, since the distance in km is divided by the travel speed (in
     * km/minute).
     * 
     * Currently, we have that the travel speed is 80km/h, i.e., 80/60 km/min.
     * @see HaversineDistance()
     */
    public static void computeDistances(double[][] timeMatrix)
    {
        HaversineDistance hd = new HaversineDistance();
        double distance = 0.0;
        for (int j = 0; j < nSites; j++)
        {
            timeMatrix[j][j] = 0.0; // diagonal
            for (int k = j+1; k < nSites; k++)
            {
                distance = hd.dist(lat[j], lon[j], lat[k], lon[k]);
                timeMatrix[j][k] = timeMatrix[k][j] = distance/_SPEED;
            }
        }
    }

    /// Greedy construction of a solution via insertion.
    /**
     *
     * Starting from a solution of type \f$\textbf{x} = ( 0, 0) \f$, where
     * \f$\textbf{x}\f$ indicates the home (starting and ending point), we
     * progressively add, in a greedy fashion, the most promising site. The
     * attractiveness of a site is measured using the followinwing score:
     * \f[\sigma_{ij} = \frac{\Delta z_i}{\Delta t_{ij}} \f]
     * where \f$\Delta z_i\f$ is the increase in the objective function value
     * associated to the inclusion of the \f$i^{th}\f$ site, and \f$\Delta
     * t_{ij}\f$ is the increase in time due to the introduction of site
     * \f$i\f$ after site \f$j\f$ (which is obviously in the current solution.)
     *
     * Therefore, for each cycle of the greedy heuristic, we determine (i)
     * which candidate should be added to the current solution (if any) and
     * (ii) in which position w.r.t. the sites in the current solution.
     *
     * The insertion mechanism is repeated as long as the budget constraint is
     * not violated. The greedy mechanism stops when either no item currently
     * not in the solution fits within the budget limits, or no items are left
     * in the `outSet`.
     */
    public static void insertionHeuristic(Solution route)
    {
        // verify the correctness of the updates
        assert ((Math.abs(route.computeTime() - route.time) < _EPSI) &&
               (Math.abs(route.computeZ() - route.Z) < _EPSI)); 

        double totTime = 0.0;
        int    totZ    = 0;

        int el = -1, pEl = -1, nEl = -1, deltaZ = 0;
        double deltaTime = 0.0, score = 0.0;
        boolean newCountry = false;

        boolean found = true;
        while (found && route.outSet.size() > 0)
        {
            Move bestMove  = new Move();
            bestMove.score = 0.0;

            found = false;
            for (int j = 0; j < route.outSet.size(); j++)
            {
                el = route.outSet.get(j); // candidate to enter inSet
                // look for best insertion in inSet
                for (int k = 0; k < route.inSet.size()-1; k++)
                {
                    pEl = route.inSet.get(k); // previous element
                    nEl = route.inSet.get(k+1); // next element 
                    deltaTime = _VISIT -timeMatrix[pEl][nEl] + timeMatrix[pEl][el] + timeMatrix[el][nEl];
                    // if item fits in budget, compute score
                    if (route.time + deltaTime <= _BUDGET)
                    {
                        deltaZ = 1 + 3*danger[el];
                        newCountry = (route.countrySet.getOrDefault(country[el], 0) == 0);
                        if (newCountry) deltaZ += 2;
                        score = (double)deltaZ/deltaTime;

                        if (score > bestMove.score)
                        {
                            found = true;
                            bestMove.update(score, deltaZ, j, k, deltaTime);
                        }
                    }
                }
            }
            if (found)
                route.insertNew(bestMove);
        }
    }
    
    /// 2-Opt move: One site from outSet moves in and one site from inSet is out.
    /** An Iterative Local Search (ILS) approach is used here. We start from an
     * initial solution, typically built via insertionHeuristic(). Next, we
     * define a <b>neighborhood</b> as the set of solutions that can be reached
     * from the incumbent swapping one site in inSet with one site from outSet.
     * In this Neighborhood, we select the solution, among the feasible ones,
     * that leads to the maximum increase in the objective function value.
     *
     * At each step of the ILS, two cases can arise:
     * <ul>
     * <li> An improving solution has been found: In this case, we set this
     * solution as incumbent and we draw a neighborhood around the new
     * incumbent. The local search restarts for a new iteration.
     * <li> An improving solution has not been found in the current
     * Neighborhood. We could apply a, e.g., Variable Neighborhood Search, to
     * increase the neighborhood size around the incumbent and restart with the
     * ILS. We did not implement such as scheme. In this version of the code,
     * when an improving solution has not been found, we stop and terminate the
     * 2-Opt heuristic.
     * </ul>
     * A classical 2-Opt move is implemented as follows. Assume a solution of the
     * following type is given:
     * \f[ \textbf{x}  = (\ldots, p_j, j, n_j, \ldots)   \f]
     * where \f$j\f$ has been chosen for exclusion from \f$\textbf{x}\f$, and
     * \f$p_j\f$ and \f$n_j\f$ are the elements preceding and followining
     * \f$j\f$, respectively. In addition, assume that element \f$k\f$ is
     * selected for insertion.
     *
     * A move that swaps \f$j\f$ with \f$k\f$ has the following
     * \f$\Delta_{kj}\f$ and the following \f$\Delta t_{kj}\f$:
     */
    public static boolean twoOptHeuristic(Solution route)
    {
        // System.out.println("2-OPT Starting Solution :: " + route.inSet);
        // System.out.println("2-OPT Starting Countries :: " + route.countrySet);
        Move bestMove = new Move();
        bestMove.deltaZ = 0;
        int elOut = -1, pEl = -1, nEl = -1, deltaZ = 0, elIn = -1, plusZ = -1;
        double deltaTime = 0.0;
        boolean newCountry = false;

        for (int j = 1; j < route.inSet.size()-1; j++)
        {
            elOut = route.inSet.get(j);
            pEl   = route.inSet.get(j-1);
            nEl   = route.inSet.get(j+1);

            // compute change in obj function value due to removal of j
            deltaZ = -(1 + 3*danger[elOut]);
            if (route.countrySet.getOrDefault(country[elOut],0) == 1)
                deltaZ -= 2; // we lose one country

            // find best candidate for inclusion
            for (int k = 0; k < route.outSet.size(); k++)
            {
                elIn = route.outSet.get(k);
                deltaTime = -timeMatrix[pEl][elOut] - timeMatrix[elOut][nEl] 
                            +timeMatrix[pEl][elIn] + timeMatrix[elIn][nEl];

                if (route.time + deltaTime > _BUDGET)
                    continue;

                plusZ = deltaZ + (1 + 3*danger[elIn]);
                // Note: two cases here: either (i) it is a new country or (ii) it is 
                // the same country of removing element and cardinality is
                // equal to 1 (we removed 2 above, we need to reintroduce 2
                // here)
                newCountry = ((route.countrySet.getOrDefault(country[elIn], 0) == 0) ||
                              ((route.countrySet.getOrDefault(country[elOut], 0)==1) && 
                              (country[elIn].equals(country[elOut]))));

                if (newCountry) // if this candidate adds a new country, count it
                    plusZ += 2;

                if (plusZ > bestMove.deltaZ) // update record of best move
                    bestMove.update(k, j, plusZ, deltaTime);
            }
        }

        if (bestMove.deltaZ > 0)
        {
            route.updateTwoOpt(bestMove);
            return false; // do not stop, local search step
        }
        else 
            return true; // stop, enhance neighborhood (use Variable Neighborhood Search)

    }


    /// Repair current solution to balance number of C and N.
    /**
     * The <i>balance</i> constraint is addressed here. The 2-Opt and the
     * construction heuristic work on a 'relaxed' version of the original
     * problem. The constraint that states that an equal number of sites of
     * type C and N should be visited (where sites of type C/N can be accounted
     * for both ways) has not been considered. Therefore, after any of the two
     * heuristic mechanisms, it is quite possible that the solution is not
     * feasible with respect to the balance constraint. 
     *
     * A repair heuristic is implemented here in two steps: We first identify
     * the imbalance. For example, as in the case of the unesco dataset, assume
     * that we have more C than N. We need a repair mechanism when:
     * \f[ |n_C - n_N | > n_{MIX} + 1   \f]
     * In other words, we allow for a maximal discrepancy of 1 site, as opposed
     * to imposing that exactly the same number of sites of the two types are
     * visited. This latter interpretation of the requirement would force the
     * solution to contain an even number of sites (which seems to be too
     * restrictive.) Here, instead, we admit a maximum difference of one site
     * for one of the two categiries. If this interpretation is not acceptable,
     * it suffices to remove one extra site from the solution produced by this
     * function.
     *
     * <ul>
     * <li> Attempt to swap one element from cateogry C currently in inSet with
     * one element of category N from outSet. We select the one with the best
     * effect on the objective function value (even if negative).
     * <li> If the repair mechanism fails, randomly remove one element from
     * inSet among those sites with Category C.
     * </ul>
     */
    public static boolean repairHeuristic(Solution route)
    {

        if (Math.abs(route.nC - route.nN) <= route.nMix +1)
            return true; // no need to repair the current solution
        else 
        {
            Move bestMove = new Move();
            int elIn = -1, elOut = -1, pEl = -1, nEl = -1, posNext = -1,  plusZ = -1, deltaZ = -1; 
            boolean newCountry;
            double swapTime, deltaTime;
            char outCat, inCat;
            bestMove.deltaZ = -_INFTY; // NOTE: this mechanism can worsen Z

            // identify current inbalance
            if (route.nC > route.nN)
            {
                outCat = 'C';
                inCat  = 'N';
            }
            else 
            {
                outCat = 'N';
                inCat  = 'C';
            }

            // identify element to be removed 
            for (int j = 1; j < route.inSet.size()-1; j++)
            {
                elOut = route.inSet.get(j);
                if (category[elOut].charAt(0) == outCat)
                {
                    pEl = route.inSet.get(j-1); 
                    nEl = route.inSet.get(j+1);
                    deltaTime = -timeMatrix[pEl][elOut] - timeMatrix[elOut][nEl] + timeMatrix[pEl][nEl];
                    deltaZ = -(1 + 3*danger[elOut]);
                    if (route.countrySet.getOrDefault(country[elOut], 0) == 1)
                        deltaZ -= 2; // we lose one country, if elOut is the only one

                    // find entering candidate, if feasible
                    for (int k = 0; k < route.outSet.size(); k++)
                    {
                        elIn = route.outSet.get(k);
                        if (category[elIn].charAt(0) == inCat) // only if of the right category
                        {
                            // it can get inserted anywhere in the current
                            // solution (this is not a swap with elOut)
                            for (int w = 0; w < route.inSet.size()-1; w++)
                            {
                                if (w == j) continue;

                                pEl = route.inSet.get(w);
                                // if next item is elOut, move one step further
                                if (route.inSet.get(w+1) == elOut)
                                {
                                    nEl = route.inSet.get(w+2);
                                    posNext = w + 2;
                                }
                                else 
                                {
                                    nEl = route.inSet.get(w+1);
                                    posNext = w + 1;
                                }
                                swapTime = deltaTime + timeMatrix[pEl][elIn] + timeMatrix[elIn][nEl] - timeMatrix[pEl][nEl];

                                if (route.time + swapTime <= _BUDGET)
                                {
                                    plusZ = deltaZ + (1 + 3*danger[elIn]);
                                    // two cases here: (i) either it is a new country, 
                                    // or it is the same country of elOut and 
                                    // we only had one of this country 
                                    newCountry = (route.countrySet.getOrDefault(country[elIn], 0) == 0 ||
                                                  ((route.countrySet.getOrDefault(country[elIn], 0)==1) &&
                                                  (country[elIn].equals(country[elOut]))));
                                    if (newCountry) plusZ += 2;

                                    if (plusZ > bestMove.deltaZ)
                                        bestMove.update(k, w, j, posNext, plusZ, swapTime);
                                }
                            }
                        }
                    }
                }
            }


            if (bestMove.deltaZ > -_INFTY) // if we found at least a feasible swap
            {
                // store here actual sited (to be used in adjusting countries)
                elOut = route.inSet.get(bestMove.Out);
                elIn  = route.outSet.get(bestMove.In);

                route.updateBalance(bestMove);
                route.Z += bestMove.deltaZ;
                // update balance counters
                if (outCat == 'C')
                {
                    route.nC -= 1;
                    route.nN += 1;
                }
                else 
                {
                    route.nN -= 1;
                    route.nC += 1;
                }
                route.time += bestMove.deltaTime;

                // update set of countries
                route.countrySet.merge(country[elOut], -1, Integer::sum);
                route.countrySet.merge(country[elIn], 1, Integer::sum);

                // verify correctness of the updates
                assert (Math.abs(route.Z - route.computeZ()) < _EPSI && 
                        Math.abs(route.time - route.computeTime()) < _EPSI);

                if (Math.abs(route.nC - route.nN) <= route.nMix + 1)
                    return true; // done, repair completed
                else
                    return false; // get out and repeat
            }
            else // we could not find a feasible swap
            {
                // activate second strategy (random removal)
                route.randomRemoveOne(outCat);
                return false;
            }
        }
    }

}

/// Manageemet of different types of moves.
/**
 * The following types of moves are managed:
 * <ul>
 * <li> Insertion for the insertionHeuristic()
 * <li> Swap for the 2-Opt heuristic
 * <li> Insertion-deletion for the repair heuristic
 * </ul>
 */
class Move
{
    int In           = -1; //!< position of element added to inSet
    int Out          = -1; //!< position of element leaving inSet
    int posIn        = -1; //!< position after which In should be inserted
    int nextIn       = -1; //!< position of next element 
    int deltaZ       = 0;  //!< change in objective function value
    double deltaTime = 0.0;//!< change in leght of tour (overall time)
    double score     =-1.0;//!< goodness score of current move

    /**
     * Update current move for insertionHeuristic().
     * One element leaves outSet and enters inSet in a specific position.
     *
     * @param score   ; ratio deltaZ/deltaTime
     * @param deltaZ  : change in obj function value
     * @param posInEl : position of element entering the solution (from outSet)
     * @param posIn   : position in which posInEl is introduced in inSet
     * @param posOut  : position of element leaving inSet
     * @param deltaTime : change in travel time 
     */
    void update(double score, int deltaZ, int posInEl, int posIn, double deltaTime)
    {
        this.score = score;
        this.deltaZ = deltaZ;
        this.In = posInEl;
        this.posIn = posIn;
        this.deltaTime = deltaTime;
    }

    void update(int inPos, int outPos, int deltaZ, double deltaTime)
    {
        this.In        = inPos;
        this.Out       = outPos;
        this.deltaZ    = deltaZ;
        this.deltaTime = deltaTime;
    }
    /**
     * Update current move for repairHeuristic.
     * One element of type outCat leaves the solution and one element of type 
     * inCat is added to the solution.
     *
     * @param posInEl : position of element entering the solution (from outSet)
     * @param posIn   : position in which posInEl is introduced in inSet
     * @param posOut  : position of element leaving inSet
     * @param posNext : position of element following posInEl
     * @param plusZ   : change in obj function value
     * @param deltaTime : change in travel time 
     */
    void update(int posInEl, int posIn, int posOut, int posNext, int deltaZ, double deltaTime)
    {
        this.In        = posInEl;
        this.posIn     = posIn;
        this.Out       = posOut;
        this.nextIn    = posNext;
        this.deltaZ    = deltaZ;
        this.deltaTime = deltaTime;
    }
}

/// Management of current solution (a route and its characteristics)
/**
 * This class is used to keep track of the characteristics of a solution, e.g,
 * the set of sites in the our of the current solution, along with the
 * cost, the overall time, the number of sites of each characteristic (C, N,
 * C/N), and the set of countries visited along the route.
 * */
class Solution
{
    ArrayList<Integer> inSet;
    ArrayList<Integer> outSet;
    int Z = 0;
    double time = 0.0;
    Map<String, Integer> countrySet;
    int nDanger = 0;
    int nC = 0;
    int nN = 0;
    int nMix = 0;
    int nCountries = 0;
    double totKm = 0.0;


    /// Initialization of a solution
    /**
     * We set the first and last sites of the route as Home and we define the
     * size of the non-dynamic data structure.
     * @param home : Index of the home site (typically 0)
     * @param nSites : Number of available sites in the dataset
     */
    Solution(int home, int nSites)
    {
        inSet = new ArrayList<Integer>(); // contains the list of sites

        // first and last elements are HOME
        inSet.add(home);
        inSet.add(home);
        outSet = new ArrayList<Integer>(); // currently, all here
        for (int j = 0; j < nSites; j++)
            if (j != home)
                outSet.add(j);

        countrySet = new HashMap<>(); // set of countries visited (empty)
    }

    /// Print solution to screen, along with its characteristics.
    void printSol()
    {
        System.out.println("---------------------------------------");
        System.out.println("              SUMMARY  ");
        System.out.println("---------------------------------------");

        System.out.println("Total Score       = " + Z);
        System.out.println("Total Sites       = " + (inSet.size()-2));
        System.out.printf( "Total Time (mins) = %.0f (vs Budget = %.0f)\n",time, Orienteering._BUDGET);
        System.out.printf( "Total Km          = %.0f\n",totKm);
        System.out.println("Nr. Endangered    = " + nDanger);
        System.out.println("Nr. Countries     = " + nCountries);
        System.out.println("Categories :: ");
        System.out.println("\t C   = "+ nC);
        System.out.println("\t N   = "+ nN);
        System.out.println("\t C/N = "+ nMix);
    }

    /// Compute full scores along all the dimensions.
    /**
     *  This is needed only at the end of the search. To foster efficiency, 
     *  we avoid recomputing solutions as much as possible. We use an update 
     *  mechanism, much faster than scanning through the entire solution.
    */
    void computeSolScores()
    {
        int el     = -1;
        int nextEl = -1;
        Z = 0; nC = 0; nN = 0; nMix = 0; nCountries = 0; nDanger = 0;
        countrySet.clear();

        // compute time (pay all visits, except for home)
        time = Orienteering._VISIT*(double)(inSet.size()-2); // exclude home twice
        // add first leg, from home to first site (cycle below starts from 1)
        time += Orienteering.timeMatrix[inSet.get(0)][inSet.get(1)];

        for (int j = 1; j < inSet.size()-1; j++)
        {
            el     = inSet.get(j);
            nextEl = inSet.get(j+1);
            time  += Orienteering.timeMatrix[el][nextEl];

            countrySet.merge(Orienteering.country[el], 1, Integer::sum);
            Z++; // +1 for a new site
            if (Orienteering.danger[el] == 1) nDanger++;
            if (Orienteering.category[el].charAt(0)  == 'C')
                nC++;
            else if (Orienteering.category[el].charAt(0)  == 'N')
                nN++;
            else 
                nMix++;
        }
        // add endangered and nr countries to the obj function
        Z += 2*countrySet.size() + 3*nDanger;
        nCountries = countrySet.size();
    }
    
    /// Compute objective function value only.
    /** This might seem redundant, but it allows to get the value of the
     * objective function without recomputing other elements that are not
     * needed.
     */
    int computeZ()
    {
        int el   = -1;
        int zVal = 0;
        nDanger  = 0;
        final Map<String, Integer> counts = new HashMap<>();
        for (int j = 1; j < inSet.size()-1; j++)
        {
            el = inSet.get(j);
            counts.merge(Orienteering.country[el], 1, Integer::sum);
            zVal++;
            if (Orienteering.danger[el] == 1) nDanger++;

        }
        zVal += 2*counts.size() + 3*nDanger;
        nCountries = counts.size();

        return zVal;
    }

    /// Computing overall time of a route  
    /** The overall time is obtained as the traveling time from home, through
     * all the sites, and back to home.
     */
    double computeTime()
    {
        int el      = -1;
        int nextEl  = -1;
        double time = Orienteering._VISIT*(double)(inSet.size()-2); // exclude home twice
        for (int j = 0; j < inSet.size()-1; j++)
        {
            el     = inSet.get(j);
            nextEl = inSet.get(j+1);
            time  += Orienteering.timeMatrix[el][nextEl];
        }
        return time;
    }
    

    /// Execute the move corresponding to the insertionHeuristic().
    /** One element is added from outSet to inSet.
     */
    void insertNew(Move move)
    {
        int pEl, nEl, el;
        pEl = inSet.get(move.posIn);
        nEl = inSet.get(move.posIn+1);
        el  = outSet.get(move.In);

        if (Orienteering.category[el].charAt(0) == 'C')
            nC += 1;
        else if (Orienteering.category[el].charAt(0) == 'N')
            nN += 1;
        else 
            nMix += 1;

        // insert element from outSet to inSet
        inSet.add(move.posIn+1, el);
        outSet.remove(move.In);
        countrySet.merge(Orienteering.country[el], 1, Integer::sum);

        time += move.deltaTime;
        Z    += move.deltaZ;

        // verify correctness of the updates
        assert ((Math.abs(computeZ() - Z) < Orienteering._EPSI) &&
                (Math.abs(computeTime() - time) < Orienteering._EPSI));
        
    }

    /// Carry out a 2-Opt swap
    /** A 2-Opt swap is carried out as follows:
     * <ul>
     * <li> insert elIn into inSet
     * <li> insert elOut into outSet
     * <li> update travel times and Objective function value
     * <li> update set of visited countries and other counters
     * </ul>
     */
    void updateTwoOpt(Move move)
    {

        int elIn = -1, elOut = -1;
        elOut = inSet.get(move.Out);
        elIn  = outSet.get(move.In);

        // insert elIn into inSet and elOut into outSet (and remove)
        inSet.add(move.Out+1, elIn);
        inSet.remove(move.Out);
        outSet.remove(move.In);
        outSet.add(elOut);

        // update countries
        countrySet.merge(Orienteering.country[elIn], 1, Integer::sum);
        countrySet.merge(Orienteering.country[elOut], -1, Integer::sum);

        // update Z and time
        Z += move.deltaZ;
        time += move.deltaTime;

        // verify correctness of the updates
        assert ((Math.abs(computeZ() - Z) < Orienteering._EPSI) &&
                (Math.abs(computeTime() - time) < Orienteering._EPSI));

    }

    /// Carry out move from repair heuristic
    /** A repair move is carried out as follows:
     * <ul>
     * <li> insert elIn into inSet in the right position, i.e., not necessarily
     * the position of outSet. This is not a swap.
     * <li> insert elOut into outSet (append it)
     * <li> update travel times and Objective function value
     * <li> update set of visited countries and other counters
     * </ul>
     */
    void updateBalance(Move move)
    {
        int pEl, nEl, elOut, elIn;

        pEl = inSet.get(move.posIn);
        nEl = inSet.get(move.nextIn);
        elOut = inSet.get(move.Out);
        elIn = outSet.get(move.In);

        // insert elIn in inSet and elOut in outSet (and remove)
        inSet.add(move.posIn+1, elIn);  // <-------- .nextIn or .posIn+1 ???
        if (move.posIn+1 <= move.Out) move.Out++; // shift by one if needed
        inSet.remove(move.Out);
        outSet.remove(move.In);
        outSet.add(elOut);
    }



    /// Randomly remove one site from a specific Category ('C' or 'N')
    /**  Generate a random number in \f$ \{1,\ldots, n-1\}\f$ (the first and
     * last are Home) and remove element associated to that position. Carry out
     * updates of Z, travel time, and other paramters.
     */
    void randomRemoveOne(char outCat)
    {
        boolean found = false;
        int pos = -1, el = -1, pEl = -1, nEl = -1;

        while (!found)
        {
            pos = 1 + Orienteering.random.nextInt(inSet.size()-2);
            el = inSet.get(pos);
            if (Orienteering.category[el].charAt(0) == outCat)
                found = true;
        }

        pEl = inSet.get(pos-1);
        nEl = inSet.get(pos+1);

        // remove el from inSet and move it to outSet
        inSet.remove(pos);
        outSet.add(el);

        // update time
        time += -Orienteering.timeMatrix[pEl][el] - Orienteering.timeMatrix[el][nEl] + Orienteering.timeMatrix[pEl][nEl];
        time -= Orienteering._VISIT;

        // update categories counters
        if (outCat == 'C')
            nC -= 1;
        else 
            nN -= 1;

        Z -= (1 + 3*Orienteering.danger[el]);
        if (countrySet.getOrDefault(Orienteering.country[el], 0) == 1)
            Z -= 2; // we lose one country
        countrySet.merge(Orienteering.country[el], -1, Integer::sum);

        // verify correctness of the updates
        assert (Math.abs(Z - computeZ()) < Orienteering._EPSI && 
        Math.abs(time - computeTime()) < Orienteering._EPSI);
    }
}

/// Store information for the instance.
/** A record is all the available information w.r.t. each site and is obtained
 * from a csv file, called `unesco.csv`.
 */
class Record
{
    int unique_number = -1; //!< unique number of a site
    int id_no         = -1; //!< id
    int danger        = -1; //!< whether a site is endangered (1) or not (0)
    double longitude  = -999; //!< longitude of the site
    double latitude   = -999; //!< latitude of the site
    String name_en; //!< name in english
    String category_short;  //!< C, N, or C/N
    String states_name_en; //!< country in english
    String iso_code; //!< country iso code

    /// Print full information of a given site
    void printRecord()
    {
        System.out.println("ID \t \t = " + id_no);
        System.out.println("Name \t \t= " + name_en);
        System.out.println("Longitude \t = " + longitude);
        System.out.println("Latitude \t = " + latitude);
        System.out.println("Category \t = " + category_short);
        System.out.println("ISO Code \t = " + iso_code);

    }
}

