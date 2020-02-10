/// Command line parser
/**
 * A number of flags can be read to define some parameters of the algorithm.
 * The command line workds as follows:
 * ~~~
 * java Orienteering -l ... -L ... -s ...
 * ~~~
 * where:
 * <ul>
 * <li> `-l` : latitude of home point. It must be a real number between -90 and
 * +90 (default : 51.533241)
 * <li> `-L` : longitude of home point. It must be a real number beween -180
 * and +180 (default : -0.104839)
 * <li> `-s` : stochastic (`-s 1`) or deterministic (`-s 00) version of the 
 * algorithm (the deterministic version fixes the seed and, therefore, always 
 * reaches the same final solution for a given set of values (longitude, latitude).
 * </ul>
 *
 * We ensure safety checks on the values of the different flags.
 */
class ParseCmdLine {
    public static void usage() 
    {
      System.out.println("Command Line Usage:   Orienteering <options>");
      System.out.println("options:       -l   latitude   [-90, 90]");
      System.out.println("options:       -L   longitude [-180, 180]");
      System.out.println("options:       -s   stochastic vs deterministic {1, 0}");
      System.out.println("Example: java Orienteering -l <latitude> -L <longitude> -s <stochastic>");
    }

    public static void readCommandLine(String[] args) 
        throws IllegalArgumentException
    {
        int i = 0, j;
        String arg;
        char flag;
        int status = 1;
        int defLongitude = 0;
        int defLatitude  = 0;
        int defStochastic = 0;

        while (i < args.length && args[i].startsWith("-")) 
        {
            arg = args[i++];

                for (j = 1; j < arg.length(); j++) {
                    flag = arg.charAt(j);
                    switch (flag) {
                    case 'l':
                        Orienteering.latitude = Double.parseDouble(args[i++]);
                        if (-90 <= Orienteering.latitude && Orienteering.latitude <= 90)
                            defLongitude = 1;
                        else 
                            defLongitude = -1;
                        break;
                    case 'L':
                        Orienteering.longitude = Double.parseDouble(args[i++]);
                        if (-180.0 <= Orienteering.longitude && Orienteering.longitude <= 180.0)
                            defLatitude = 1;
                        else
                            defLatitude = -1;
                        break;
                    case 's':
                        Orienteering.stochastic = Integer.parseInt(args[i++]);
                        if ((Orienteering.stochastic != 1) && (Orienteering.stochastic != 0))
                            defStochastic = -1;
                        break;
                    case 'h':
                        throw new IllegalArgumentException ();
                    default:
                        System.err.println("ParseCmdLine: Illegal option " + flag);
                        break;
                    }
                }
        }
        if (defLongitude == -1 || defLatitude == -1 || defStochastic == -1)
            throw new IllegalArgumentException ();
    }
}
