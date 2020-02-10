import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.Scanner;

/// Management of I/O from/to the disk
/** The input/output is managed via a single csv file, stored under the folder
 * data, and called `unesco.csv`. This class creates an instance of Record()
 * and stores the full information into a list of records, called DataFrame.
 */
public class ReadCSV
{
	//Delimiters used in the CSV file
	private static final String COMMA_DELIMITER = ",";

    // dataset for the challenge
	private static final String namefile = "data/unesco.csv";
	
    /// Read data from disk
	public static void readData()
	{
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(namefile));

			String line = reader.readLine();
			line = reader.readLine();
			while (line != null) {
                Record record = new Record();
                String[] splittedLine = line.split(COMMA_DELIMITER);
                record.unique_number = Integer.valueOf(splittedLine[0]);
                record.id_no = Integer.valueOf(splittedLine[1]);
                record.name_en = splittedLine[2];
                record.danger = Integer.valueOf(splittedLine[3]);
                record.longitude = Double.valueOf(splittedLine[4]);
                record.latitude = Double.valueOf(splittedLine[5]);
                record.category_short = splittedLine[6];
                record.states_name_en = splittedLine[7];
                record.iso_code = splittedLine[8];

                Orienteering.DataFrame.add(record);

                // read next line
				line = reader.readLine();
			}
			reader.close();


		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}



