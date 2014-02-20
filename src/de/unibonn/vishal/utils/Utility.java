/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.unibonn.vishal.utils;

import java.io.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.nio.channels.FileChannel;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JTextArea;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.JTextComponent;

/**
 *
 * @author tckb
 */
public class Utility {

    // Public Meta-Inf
    public final static String About = " Utility toolbox";
    public final static String Ver = "0.1 - beta";
    public final static String Author = "Chandra Tungathurthi";
    // Public constants
    public final static String WORD_BREAK = " ";
    public final static String LINE_BREAK = System.getProperty("line.separator");
    public final static String FILE_SEPERATOR = System.getProperty("file.separator");
    public final static String OS_NAME = System.getProperty("os.name");
    public final static String OS_ARC = System.getProperty("os.arch");
    public final static String J_CLSPTH = System.getProperty("java.class.path");
    public final static String USER_HME = System.getProperty("user.home");
    public final static String USER_DIR = System.getProperty("user.dir");
    /**
     * Available computing power in terms of available processors
     */
    public final static int compCores = Runtime.getRuntime().availableProcessors();
    public final static FileNameExtensionFilter wavFileFilter = new FileNameExtensionFilter("(*.wav) Microsoft Wave files", "wav");
    public final static FileNameExtensionFilter txtFileFilter = new FileNameExtensionFilter("(*.txt ) Text files", "txt");
    // Private stuff
    private static final Logger mylogger = Logger.getLogger("com.lia.core");

    private static void copy12(File file1, File file2) {
        try {
            mylogger.log(Level.INFO, "Copy12:[0}] -> [{1}]{ ...", new Object[]{file1.getAbsolutePath(), file2.getAbsolutePath()});

            FileChannel in = (new FileInputStream(file1)).getChannel();
            FileChannel out = (new FileOutputStream(file2)).getChannel();
            in.transferTo(0, file1.length(), out);
            in.close();
            out.close();
            mylogger.info("... done}");
        } catch (Exception ex) {
            mylogger.log(Level.SEVERE, "Error while copying: ", ex);

        }

    }

    // Static methods 
    /**
     *
     * @param fromFile
     * @param toFile
     * @return
     */
    public static File copyToFolderAs(File fromFile, File folder, String asName) {
        mylogger.log(Level.INFO, "Copying:[{0}] file to folder [{1}]  as {2}", new Object[]{fromFile.getAbsolutePath(), folder.getAbsolutePath(), asName});

        if (!folder.isDirectory() || !folder.exists()) {
            mylogger.log(Level.SEVERE, "Invalid folder {0} ! COPY FAILED", new Object[]{folder.getAbsolutePath()});

            return null;
        }
        if (!fromFile.exists()) {
            mylogger.log(Level.SEVERE, "File {0} is doesn't exists! COPY FAILED", new Object[]{fromFile.getAbsolutePath()});

            return null;
        } else {
            try {
                File newFile = new File(folder, asName);
                if (!newFile.createNewFile()) {
                    mylogger.log(Level.SEVERE, "File {0} creation failed!", new Object[]{newFile.getAbsolutePath()});

                    return null;
                } else {
                    copy12(fromFile, newFile);
                    mylogger.info("Copy completed");
                    return newFile;
                }



            } catch (IOException ex) {
                mylogger.log(Level.SEVERE, "Something went wrong; error while copying: ", ex);
            }
        }
        return null;

    }

    /**
     *
     * @param fileToCopy
     * @param folder
     * @return
     */
    public static boolean copyToFolder(File fileToCopy, File folder) {

        if (copyToFolderAs(fileToCopy, folder, fileToCopy.getName()) != null) {
            return true;
        } else {
            return false;
        }
    }

    public static File makeDuplicate(File thisFile) {
        mylogger.log(Level.INFO, "Duplicating {0}", new Object[]{thisFile.getAbsolutePath()});
        return copyToFolderAs(thisFile, thisFile.getParentFile(), stripExtension(thisFile) + "_copy" + getExtension(thisFile));

    }

    public static String getSafePath(File file) {
        return "\"" + file.getPath() + "\"";
    }

    public static String getSafeName(File file) {
        return file.getName().trim().replace(".", "_");
    }

    // Apparently, fast load text file ? "readin" function in original text
    // http://users.cs.cf.ac.uk/O.F.Rana/jdc/swing-nov7-01.txt
    synchronized public static void loadFileToPane(String fname, JTextComponent pane) {
        FileReader fr = null;
        try {
            fr = new FileReader(fname);
            pane.read(fr, null);
            fr.close();
        } catch (IOException ex) {
            mylogger.log(Level.SEVERE, "Error while loading: ", ex);
        } finally {
            try {
                fr.close();
            } catch (IOException ex) {
                mylogger.log(Level.SEVERE, "Error while loading: ", ex);
            }
        }

    }

    synchronized public static void saveStringToFile(String text, File fname) {
        try {
            mylogger.log(Level.INFO, "Loading: [{0}] ->file: [{1}]", new Object[]{text, fname.getAbsolutePath()});

            FileWriter fr = new FileWriter(fname);
            fr.write(text);

            fr.close();
        } catch (IOException ex) {
            mylogger.log(Level.SEVERE, "Error while loading: ", ex);
        }

    }

    synchronized public static void loadTextToPane(String text, JTextComponent pane, boolean append) {

        if (append) {
            ((JTextArea) pane).append(text);
        } else {
            ((JTextArea) pane).setText(text);
        }
    }

    public static String stripExtension(File file) {
        mylogger.log(Level.INFO, "Stripping extension: {0}", file.getName());
        int i = file.getName().lastIndexOf(".");
        String s = file.getName().substring(0, i);
        mylogger.log(Level.FINE, ". location{0}", i);
        mylogger.log(Level.FINE, "name {0}", s);
        return s;
    }

    public static String getExtension(File file) {
        String fname = file.getName();
        mylogger.log(Level.INFO, "Getting extension: {0}", fname);
        String fNoext = stripExtension(file);
        return fname.substring(fname.indexOf(fNoext) + fNoext.length());
    }

    /**
     * Returns individual words in 'f' as a Array of Strings
     *
     * @param f - file to be read
     * @return Array of Strings
     */
    public static String[] getWordsInFile(File f) {

        return readFileAsLongString(f).split(WORD_BREAK);
    }

    /**
     *
     * @param f - file name
     * @return File contents as string
     */
    private static String readFile(File f, boolean preserveLineBreaks) {

        StringBuilder content = new StringBuilder();
        String line = null;
        BufferedReader br = null;
        try {

            mylogger.log(Level.INFO, "Reading: [{0}] File: [{1}]", new Object[]{f.getName(), f.getAbsolutePath()});
            br = new BufferedReader(new FileReader(f));
            while ((line = br.readLine()) != null) {
                if (preserveLineBreaks) {
                    content.append(line).append(LINE_BREAK);

                } else {
                    content.append(line).append(WORD_BREAK);
                }
            }

        } catch (FileNotFoundException ex) {
            mylogger.log(Level.SEVERE, "Error:", ex);
        } finally {

            try {

                br.close();
            } catch (IOException ex) {
                mylogger.log(Level.SEVERE, "Error:", ex);
            }
            return content.toString();
        }
    }

    public static String readFileAsString(String fname) {


        return readFile(new File(fname), true);
    }

    public static String readFileAsString(File f) {


        return readFile(f, true);
    }

    public static String readFileAsLongString(String fname) {


        return readFile(new File(fname), false);
    }

    public static String readFileAsLongString(File f) {


        return readFile(f, false);
    }

    /**
     * Returns offset of the given word in src the component
     *
     * @param src source component
     * @param word desired word
     * @return position offset
     */
    public static int searchInTxComp(JTextComponent src, String word) {
        int firstOffset = -1;


        if (word == null || word.isEmpty()) {
            return -1;
        }

        // Look for the word we are given - insensitive searchInTxComp
        String content = null;
        try {
            Document d = src.getDocument();

            content = d.getText(0, d.getLength()).toLowerCase();
        } catch (BadLocationException e) {
            // Cannot happen
            return -1;
        }

        word = word.toLowerCase();
        int lastIndex = 0;
        int wordSize = word.length();




        while ((lastIndex = content.indexOf(word, lastIndex)) != -1) {
            int endIndex = lastIndex + wordSize;

            if (firstOffset == -1) {
                firstOffset = lastIndex;
            }
            lastIndex = endIndex;
        }

        return firstOffset;
    }

    /**
     *
     * @param transcriptFile
     * @return Hashmap of (timestamp, transcriptword) pair
     */
    public static HashMap<Double, String> getTrList(File transcriptFile) {
        HashMap<Double, String> trList = new HashMap<Double, String>();

        String line = null;
        BufferedReader br = null;
        try {

            mylogger.log(Level.INFO, "Parsing transcript file:{0}", new Object[]{transcriptFile.getName()});
            br = new BufferedReader(new FileReader(transcriptFile));
            while ((line = br.readLine()) != null) {

                String parts[] = line.split(" ");

                Double ts = Double.parseDouble(parts[0]);

                String word = parts[1];
                trList.put(ts * 1000, word); // store the timestamps as ms

                //  System.out.println("Contains key:" + trList.containsKey(ts) + ":" + trList.get(ts));
            }
            // System.out.println(trList);
            mylogger.log(Level.INFO, "Parsing done; trlist size:{0}", trList.size());

        } catch (FileNotFoundException ex) {
            mylogger.log(Level.SEVERE, "Error:", ex);
        } finally {

            try {

                br.close();
            } catch (IOException ex) {
                mylogger.log(Level.SEVERE, "Error:", ex);
            }


            return trList;
        }


    }

    public static File createTmpFile(String prx, String sfx) {
        try {
            return File.createTempFile(prx, sfx);
        } catch (IOException ex) {
            mylogger.log(Level.SEVERE, "Error: Can not create temp file ", ex);
            return null;
        }
    }

    public static File createTmpFile() {
        try {
            return File.createTempFile("Utility_temp_file-" + String.valueOf(Utility.tic()), ".tmp");
        } catch (IOException ex) {
            mylogger.log(Level.SEVERE, "Error: Can not create temp file ", ex);
            return null;
        }
    }

    public static double adjDecimalSep(double value) {
        mylogger.log(Level.FINE, "Adjusting valu {0}", new Object[]{value});
        // Bug fix: 
//        For german locale 3.333 => 3,33 which would raise an error 
        DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.getDefault());
        otherSymbols.setDecimalSeparator('.');
        otherSymbols.setGroupingSeparator(',');
        // Bug fix
        DecimalFormat newFormat = new DecimalFormat();
        newFormat.setDecimalFormatSymbols(otherSymbols);
        return Double.valueOf(newFormat.format(value));

    }

    public static double roundDouble(Double val, int precision) {

        Double dVal = adjDecimalSep(val);

        BigDecimal bigNumber1 = new BigDecimal(dVal, MathContext.UNLIMITED); // dec 64 - > double precision
        return bigNumber1.setScale(precision, BigDecimal.ROUND_CEILING).doubleValue();

    }

    public static long tic() {
        return System.nanoTime();
    }

    public static double toc(long tic) {
        return (double) (tic() - tic) / Math.pow(10, 9);
    }

    public static void saveObjectToFile(File thatFile, Object... thisObject) {
        mylogger.log(Level.INFO, "Saving object{0} to file: {1}", new Object[]{thisObject.hashCode(), thatFile.getAbsolutePath()});
        ObjectOutputStream stream;
        try {
            thatFile.delete();
            if(thatFile.createNewFile()){
            stream = new ObjectOutputStream(new FileOutputStream(thatFile));
            stream.writeObject(thisObject);
            stream.close();
            }else{
                mylogger.warning("File creation failed! try again");
            }
                
        } catch (IOException ex) {
            mylogger.log(Level.SEVERE, "Error: Can not save the object!", ex.getMessage());
        }

    }

    public static Object[] getObjectFromFile(File thatFile) {
        mylogger.log(Level.INFO, "Extracting object from file: {0}", thatFile.getAbsolutePath());

        Object[] thatObject = null;
        try {
            ObjectInputStream stream = new ObjectInputStream(new FileInputStream(thatFile));
            thatObject = (Object[]) stream.readObject();
            stream.close();
        } catch (Exception ex) {
            mylogger.log(Level.SEVERE, "Error: Can not extract the object from file", ex.getMessage());
        }
        return thatObject;
    }

    private Utility() {
    }

    // Shamefully borrowed from Stackoverflow!
    // http://stackoverflow.com/questions/109383/how-to-sort-a-mapkey-value-on-the-values-in-java
    /**
     * @author Maxim Veksler
     */
    public static class MapUtils {

        public static int countUniqueValues(Map<String, Integer> map) {
            int cnt = 0;
            ArrayList<Integer> valueList = new ArrayList<Integer>(map.values());
            for (Integer val : valueList) {
                if (val != -Integer.MAX_VALUE) {
                    Collections.replaceAll(valueList, val, -Integer.MAX_VALUE);
                    cnt++;
                }
            }

            return cnt;
        }

        /**
         * Sort a map by it's keys in ascending order.
         *
         * @return new instance of {@link LinkedHashMap} contained sorted
         * entries of supplied map.
         * @author Maxim Veksler
         */
        public static <K, V> LinkedHashMap<K, V> sortMapByKey(final Map<K, V> map) {
            return sortMapByKey(map, SortingOrder.ASCENDING);
        }

        /**
         * Sort a map by it's values in ascending order.
         *
         * @return new instance of {@link LinkedHashMap} contained sorted
         * entries of supplied map.
         * @author Maxim Veksler
         */
        public static <K, V> LinkedHashMap<K, V> sortMapByValue(final Map<K, V> map) {
            return sortMapByValue(map, SortingOrder.ASCENDING);
        }

        /**
         * Sort a map by it's keys.
         *
         * @param sortingOrder {@link SortingOrder} enum specifying requested
         * sorting order.
         * @return new instance of {@link LinkedHashMap} contained sorted
         * entries of supplied map.
         * @author Maxim Veksler
         */
        public static <K, V> LinkedHashMap<K, V> sortMapByKey(final Map<K, V> map, final SortingOrder sortingOrder) {
            Comparator<Map.Entry<K, V>> comparator = new Comparator<Entry<K, V>>() {
                public int compare(Entry<K, V> o1, Entry<K, V> o2) {
                    return comparableCompare(o1.getKey(), o2.getKey(), sortingOrder);
                }
            };

            return sortMap(map, comparator);
        }

        /**
         * Sort a map by it's values.
         *
         * @param sortingOrder {@link SortingOrder} enum specifying requested
         * sorting order.
         * @return new instance of {@link LinkedHashMap} contained sorted
         * entries of supplied map.
         * @author Maxim Veksler
         */
        public static <K, V> LinkedHashMap<K, V> sortMapByValue(final Map<K, V> map, final SortingOrder sortingOrder) {
            Comparator<Map.Entry<K, V>> comparator = new Comparator<Entry<K, V>>() {
                public int compare(Entry<K, V> o1, Entry<K, V> o2) {
                    return comparableCompare(o1.getValue(), o2.getValue(), sortingOrder);
                }
            };

            return sortMap(map, comparator);
        }

        @SuppressWarnings("unchecked")
        private static <T> int comparableCompare(T o1, T o2, SortingOrder sortingOrder) {
            int compare = ((Comparable<T>) o1).compareTo(o2);

            switch (sortingOrder) {
                case ASCENDING:
                    return compare;
                case DESCENDING:
                    return (-1) * compare;
            }

            return 0;
        }

        /**
         * Sort a map by supplied comparator logic.
         *
         * @return new instance of {@link LinkedHashMap} contained sorted
         * entries of supplied map.
         * @author Maxim Veksler
         */
        public static <K, V> LinkedHashMap<K, V> sortMap(final Map<K, V> map, final Comparator<Map.Entry<K, V>> comparator) {
            // Convert the map into a list of key,value pairs.
            List<Map.Entry<K, V>> mapEntries = new LinkedList<Map.Entry<K, V>>(map.entrySet());

            // Sort the converted list according to supplied comparator.
            Collections.sort(mapEntries, comparator);

            // Build a new ordered map, containing the same entries as the old map.  
            LinkedHashMap<K, V> result = new LinkedHashMap<K, V>(map.size() + (map.size() / 20));
            for (Map.Entry<K, V> entry : mapEntries) {
                // We iterate on the mapEntries list which is sorted by the comparator putting new entries into 
                // the targeted result which is a sorted map. 
                result.put(entry.getKey(), entry.getValue());
            }

            return result;
        }

        private MapUtils() {
        }

        /**
         * Sorting order enum, specifying request result sort behavior.
         *
         * @author Maxim Veksler
         *
         */
        public static enum SortingOrder {

            /**
             * Resulting sort will be from smaller to biggest.
             */
            ASCENDING,
            /**
             * Resulting sort will be from biggest to smallest.
             */
            DESCENDING
        }
    }

    public static class UI {

        public static File getFile(JComponent parent) {
            JFileChooser jfc = new JFileChooser();
            jfc.setDialogTitle("Utility: FileChooser");
            jfc.showOpenDialog(parent);

            return jfc.getSelectedFile();

        }

        public static File saveFile(JComponent parent) {
            File f = null;
            try {
                JFileChooser jfc = new JFileChooser();
                jfc.setDialogTitle("Utility: FileSaver");
                jfc.showSaveDialog(parent);
                f = jfc.getSelectedFile();
                f.createNewFile();
            } catch (IOException ex) {
                mylogger.log(Level.SEVERE, "Error: Can not save file: {0}", ex.getMessage());
            } finally {
                return f;

            }
        }

        public static File getFile(JComponent parent, FileNameExtensionFilter filter) {
            JFileChooser jfc = new JFileChooser();
            jfc.addChoosableFileFilter(filter);
            jfc.setDialogTitle("Utility: FileChooser");
            jfc.showOpenDialog(parent);
            return jfc.getSelectedFile();

        }

        public static void showInfoMessage(JComponent parent, String message) {
            JOptionPane.showMessageDialog(parent, message, "Utility: Info", JOptionPane.INFORMATION_MESSAGE);

        }
    }
}
// -DEAD CODE -
//    public static File createOutFile(DIR type, String parent, String child) throws IOException {
//
//        File parentDir = new File((File) curOutDir.get(type), parent);
//        if (parentDir.mkdir()) {
//            File childFile = new File(parentDir, child); // then create child inside "parent"
//            if (childFile.createNewFile()) {
//                return childFile;
//            }
//        } // create directory inside "type"
//
//        return null;
//
//    }

