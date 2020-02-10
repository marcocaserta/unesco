/// Implementation Haversine Distance Algorithm between two points.
/**
 * @see <a href="https://bigdatanerd.wordpress.com/2011/11/03/java-implementation-of-haversine-formula-for-distance-calculation-between-two-points/">See source here.</a>
 * <ul>
 * <li> R = earth’s radius (mean radius = 6,371km)
 * <li> Δlat = lat2− lat1
 * <li> Δlong = long2− long1
 * <li> a = sin²(Δlat/2) + cos(lat1).cos(lat2).sin²(Δlong/2)
 * <li> c = 2.atan2(√a, √(1−a))
 * <li> d = R.c
 * </ul>

 * @param lat1 : latitude point 1 
 * @param lon1 : longitude point 1 
 * @param lat2 : latitude point 2
 * @param lon2 : latitude point 2 
 *
 * @return distance: distance in km between the two points
 */
public class HaversineDistance {
 
    public static double dist(double lat1, double lon1, double lat2, double lon2) {
        final int R = 6371; // Radious of the earth
        Double latDistance = toRad(lat2-lat1);
        Double lonDistance = toRad(lon2-lon1);
        Double a = Math.sin(latDistance / 2) * Math.sin(latDistance / 2) + 
                   Math.cos(toRad(lat1)) * Math.cos(toRad(lat2)) * 
                   Math.sin(lonDistance / 2) * Math.sin(lonDistance / 2);
        Double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
        Double distance = R * c;

        return distance;
         
        /* System.out.println("The distance between two lat and long is::" + distance); */
 
    }
     
    private static Double toRad(Double value) {
        return value * Math.PI / 180;
    }
 
}
