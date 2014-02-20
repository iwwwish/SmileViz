/*
 * Copyright (C) 2014. EMBL, European Bioinformatics Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.unibonn.vishal.utils;

import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.vecmath.Point2d;
import net.sf.jniinchi.INCHI_RET;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.exception.NoSuchAtomTypeException;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fragment.IFragmenter;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.CMLWriter;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.atomic.AtomDegreeDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.AtomHybridizationDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.AtomValenceDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.CovalentRadiusDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.VdWRadiusDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

/**
 *
 * @author Vishal Siramshetty <vishal[at]ebi.ac.uk>
 */
public class CDKUtil {

    /**
     * Provides methods that can be implemented on CDK Atom.
     */
    public static class Atoms {

        /**
         * Returns the non-hydrogen neighbour atom count for a given atom in a
         * molecule.
         *
         * @param atom
         * @param mol
         * @return
         */
        public static int getDegree(IAtom atom, IAtomContainer mol) {

            DescriptorValue result = new AtomDegreeDescriptor().calculate(atom, mol);
            return Integer.parseInt(result.getValue().toString());
        }

        /**
         * Returns the hybridization state for a given atom in a molecule.
         *
         * @param atom
         * @param mol
         * @return
         */
        public static String getHybridization(IAtom atom, IAtomContainer mol) {

            DescriptorValue result = new AtomHybridizationDescriptor().calculate(atom, mol);

            return result.getValue().toString();

        }

        /**
         * Returns the valency for a given atom in a molecule.
         *
         * @param atom
         * @param mol
         * @return
         */
        public static int getValency(IAtom atom, IAtomContainer mol) {

            DescriptorValue result = new AtomValenceDescriptor().calculate(atom, mol);
            return Integer.parseInt(result.getValue().toString());
        }

        /**
         * Returns the covalent radius for a given atom in a molecule.
         *
         * @param atom
         * @param mol
         * @return
         */
        public static Double getCovalentRadius(IAtom atom, IAtomContainer mol) {
            try {
                DescriptorValue result = new CovalentRadiusDescriptor().calculate(atom, mol);
                return Double.parseDouble(result.getValue().toString());
            } catch (IOException | ClassNotFoundException ex) {
                Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
            }
            return null;
        }

        /**
         * Returns the Van der Waals radius for a given atom in a molecule.
         *
         * @param atom
         * @param mol
         * @return
         */
        public static Double getVanderWaalsRadius(IAtom atom, IAtomContainer mol) {
            try {
                DescriptorValue result = new VdWRadiusDescriptor().calculate(atom, mol);
                return Double.parseDouble(result.getValue().toString());
            } catch (IOException | ClassNotFoundException ex) {
                Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
            }
            return null;
        }

        private Atoms() {

        }
    }

    /**
     * Provides methods that can be implemented on a CDK AtomContainer.
     */
    public static class Molecule {

        /**
         * Returns string representation of the molecular formula for a given
         * molecule.
         *
         * @param mol
         * @return
         */
        public static String getMolecularFormula(IAtomContainer mol) {

            IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(mol);

            return MolecularFormulaManipulator.getString(formula);
        }

        /**
         * Returns XlogP for a given molecule.
         *
         * @param mol
         * @return
         */
        public static Double getXLogP(IAtomContainer mol) {

            DescriptorValue result = new XLogPDescriptor().calculate(mol);

            return Math.round(Double.parseDouble(result.getValue().toString()) * 100.0) / 100.0;

        }

        /**
         * Returns the topological surface area for a given molecule.
         *
         * @param mol
         * @return
         */
        public static Double getTPSA(IAtomContainer mol) {

            DescriptorValue result = new TPSADescriptor().calculate(mol);

            return Math.round(Double.parseDouble(result.getValue().toString()) * 100.0) / 100.0;
        }

        /**
         * Returns the exact molecular weight for a given molecule.
         *
         * @param mol
         * @return
         */
        public static Double getExactMass(IAtomContainer mol) {

            int hCount = AtomContainerManipulator.getImplicitHydrogenCount(mol);

            double avgMass = AtomContainerManipulator.getNaturalExactMass(mol);
            avgMass = Math.round((avgMass) * 100.0) / 100.0;

            return avgMass + hCount;
        }

        /**
         * Returns SMILES for a given molecule.
         *
         * @param mol
         * @return
         * @throws CDKException
         */
        public static String getSMILES(IAtomContainer mol) throws CDKException {
            return new SmilesGenerator().create(mol);
        }

        /**
         * Returns the unique SMILES representation for given molecule.
         *
         * @param mol
         * @return
         * @throws InvalidSmilesException
         * @throws CDKException
         */
        public static String getUniqueSMILES(IAtomContainer mol) throws InvalidSmilesException, CDKException {

            String smile = SmilesGenerator.unique().create(mol);

            return smile;
        }

        public static String getBemisMurckoSMILES(IAtomContainer mol) throws CDKException {
            IFragmenter fg = new MurckoFragmenter();
            fg.generateFragments(mol);
            String[] bmScaffs = fg.getFragments();
            return getLongestString(bmScaffs);
        }

        private static String getLongestString(String[] fragments) {
            String bigstring = null;
            int maxlength = 0;
            for (String max : fragments) {
                if (maxlength < max.length()) {
                    maxlength = max.length();
                    bigstring = max;
                }
            }
            return bigstring;
        }

        /**
         * Returns the absolute SMILES representation for given molecule.
         *
         * @param mol
         * @return
         * @throws InvalidSmilesException
         * @throws CDKException
         */
        public static String getAbsoluteSMILES(IAtomContainer mol) throws InvalidSmilesException, CDKException {

            String smile = SmilesGenerator.absolute().create(mol);

            return smile;
        }

        /**
         * Returns the aromatic SMILES representation for given molecule.
         *
         * @param mol
         * @return
         * @throws org.openscience.cdk.exception.CDKException
         */
        public static String getAromaticSMILES(IAtomContainer mol) throws CDKException {
            String smile = SmilesGenerator.unique().aromatic().create(mol);
            return smile;
        }

        /**
         * Returns InChI for a given molecule.
         *
         * @param mol
         * @return
         */
        public static String getInChI(IAtomContainer mol) {
            try {
                InChIGenerator generator = InChIGeneratorFactory.getInstance().getInChIGenerator(mol);
                if (generator.getReturnStatus() == INCHI_RET.OKAY) {
                    return generator.getInchi();
                } else if (generator.getMessage().startsWith("Accepted")) {
                    AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
                    return generator.getInchi();
                } else {
                    System.out.println("InChI was not generated due to this error: " + generator.getMessage());
                }

            } catch (CDKException ex) {
                Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
            }
            return null;
        }

        /**
         * Returns InChI Key for a given molecule.
         *
         * @param mol
         * @return
         */
        public static String getInChIKey(IAtomContainer mol) {
            try {
                InChIGenerator generator = InChIGeneratorFactory.getInstance().getInChIGenerator(mol);
                return generator.getInchiKey();
            } catch (CDKException ex) {
                Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
            }
            return null;
        }

        /**
         * Returns the rotatable bond count for a given molecule.
         *
         * @param mol
         * @return
         */
        public static int getRotatableBondCount(IAtomContainer mol) {

            DescriptorValue result = new RotatableBondsCountDescriptor().calculate(mol);

            return Integer.parseInt(result.getValue().toString());

        }

        /**
         * Returns number of aromatic atoms in a given molecule.
         *
         * @param mol
         * @return
         */
        public static int getAromaticAtomCount(IAtomContainer mol) {

            DescriptorValue result = new AromaticAtomsCountDescriptor().calculate(mol);

            return Integer.parseInt(result.getValue().toString());

        }

        /**
         * Returns number of aromatic bonds in a given molecule.
         *
         * @param mol
         * @return
         */
        public static int getAromaticBondCount(IAtomContainer mol) {

            DescriptorValue result = new AromaticBondsCountDescriptor().calculate(mol);

            return Integer.parseInt(result.getValue().toString());
        }

        /**
         * Returns the number of hydrogen bond acceptors in a given molecule.
         *
         * @param mol
         * @return
         */
        public static int getHBondAcceptorCount(IAtomContainer mol) {

            DescriptorValue result = new HBondAcceptorCountDescriptor().calculate(mol);

            return Integer.parseInt(result.getValue().toString());
        }

        /**
         * Returns the number of hydrogen bond donors in a given molecule.
         *
         * @param mol
         * @return
         */
        public static int getHBondDonorCount(IAtomContainer mol) {

            DescriptorValue result = new HBondDonorCountDescriptor().calculate(mol);

            return Integer.parseInt(result.getValue().toString());
        }

        /**
         * Returns the number of heteroatoms (atoms other than Carbon and
         * Hydrogen) in a given molecule.
         *
         * @param mol
         * @return
         */
        public static int getHeteroAtomCount(IAtomContainer mol) {

            int count = 0;
            for (IAtom atom : mol.atoms()) {
                if (atom.getAtomicNumber() != 1 && atom.getAtomicNumber() != 6) {
                    count++;
                }
            }

            return count;

        }

        /**
         * Returns a map containing properties for a given molecule.
         *
         * @param mol
         * @return
         * @throws CDKException
         */
        public static TreeMap<Object, Object> getProperties(IAtomContainer mol) throws CDKException {

            TreeMap<Object, Object> properties = new TreeMap<>();
            properties.put("Molecular Formula", getMolecularFormula(mol));
            properties.put("SMILES", getSMILES(mol));
            properties.put("InChI", getInChI(mol));
            properties.put("InChI Key", getInChIKey(mol));
            properties.put("Molecular Weight", getExactMass(mol));
            properties.put("XLog P", getXLogP(mol));
            properties.put("TPSA", getTPSA(mol));
            properties.put("Rotatable Bond Count", getRotatableBondCount(mol));
            properties.put("Hydrogen Bond Acceptor Count", getHBondAcceptorCount(mol));
            properties.put("Hydrogen Bond Donor Count", getHBondDonorCount(mol));
            properties.put("Aromatic Atom Count", getAromaticAtomCount(mol));
            properties.put("Aromatic Bond Count", getAromaticBondCount(mol));

            return properties;

        }

        /**
         * Returns true if a given molecule obeys Lipinski's Rule of five or
         * false if it does not obey at least one rule.
         *
         * @param mol
         * @return
         */
        public static boolean obeysRuleOfFive(IAtomContainer mol) {

            DescriptorValue result = new RuleOfFiveDescriptor().calculate(mol);
            return Integer.parseInt(result.getValue().toString()) <= 0;
        }

        /**
         * Returns the MACCS (166 bit) fingerprint for the given molecule.
         *
         * @param mol
         * @return
         * @throws CDKException
         */
        public static String getMACCSFP(IAtomContainer mol) throws CDKException {

            MACCSFingerprinter maccsFPer = new MACCSFingerprinter();
            BitSet maccsFP = maccsFPer.getBitFingerprint(mol).asBitSet();

            return maccsFP.toString();

        }

        /**
         * Returns the string representation of MACCS fingerprint for given
         * molecule.
         *
         * @param mol
         * @return {1,2,3,4,.....,166(0 or 1)}
         * @throws CDKException
         */
        public static String getMACCSFPString(IAtomContainer mol) throws CDKException {

            MACCSFingerprinter maccsFPer = new MACCSFingerprinter();
            BitSet maccsFP = maccsFPer.getBitFingerprint(mol).asBitSet();

            int[] fp = new int[166];

            for (int i = 0; i < fp.length; i++) {
                if (maccsFP.get(i) == true) {
                    fp[i] = 1;
                } else {
                    fp[i] = 0;
                }

            }

            return Arrays.toString(fp);
        }

        /**
         * Returns true if two molecules are identical by comparing their
         * InChIs.
         *
         * @param mol1
         * @param mol2
         * @return
         */
        public static boolean moleculesIdentical(IAtomContainer mol1, IAtomContainer mol2) {

            String inchi1 = getInChI(mol1);
            String inchi2 = getInChI(mol2);
            return inchi1.equalsIgnoreCase(inchi2);
        }

        /**
         * Returns Tanimoto similarity between mol1 and mol2.
         *
         * @param mol1
         * @param mol2
         * @return
         */
        public static Double getTanimotoSimilarity(IAtomContainer mol1, IAtomContainer mol2) {
            MACCSFingerprinter maccsFPer = new MACCSFingerprinter();
            int a = 0, b = 0, c = 0;
            double sim = 0;
            try {
                BitSet maccsFP1 = maccsFPer.getBitFingerprint(mol1).asBitSet();
                BitSet maccsFP2 = maccsFPer.getBitFingerprint(mol2).asBitSet();
                a = maccsFP1.cardinality();
                b = maccsFP2.cardinality();
                c = 0;
                for (int i = 0; i <= maccsFP1.length(); i++) {
                    for (int j = i; j < maccsFP2.length(); j++) {
                        if (maccsFP1.get(i) == true && maccsFP2.get(j) == true) {
                            c++;
                        }
                        i++;
                    }
                }
                sim = c / (a + b - c);
                return sim;

            } catch (CDKException ex) {
                System.err.println("CDKException: " + ex.getMessage());
            }
            return null;
        }

        /**
         * Returns the list of atoms present in the shortest path between given
         * two atoms.
         *
         * @param mol
         * @param atom1
         * @param atom2
         * @return
         */
        public static List<IAtom> getAtomsInShortestPath(IAtomContainer mol, IAtom atom1, IAtom atom2) {

            return PathTools.getShortestPath(mol, atom1, atom2);

        }

        private Molecule() {

        }
    }

    /**
     * Provides methods for storing and retrieving molecules in various IO
     * formats.
     */
    public static class IO {

        /**
         * Saves the molecule object as SDfile to a specified directory or to
         * home directory, if unspecified.
         *
         * @param molecule
         * @param path
         */
        public static void saveAsSDF(IAtomContainer molecule, String path) {

            if (path == null) {
                path = System.getProperty("user.home");
            }
            try {
                OutputStream out = new FileOutputStream(path.concat("/output.sdf"));
                SDFWriter writer = new SDFWriter(out);
                writer.write(molecule);
                writer.close();
            } catch (FileNotFoundException | CDKException ex) {
                Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                System.err.println("IOException" + ex.getMessage());
            }
        }

        /**
         * Saves the molecule object as SDfile, along with its properties to a
         * specified directory or to home directory, if unspecified.
         *
         * @param molecule
         * @param path
         */
        public static void saveAsSDFWithProperties(IAtomContainer molecule, String path) {

            try {
                if (path == null) {
                    path = System.getProperty("user.home") + File.separator + "output.sdf";
                }
                OutputStream out = new FileOutputStream(path);
                TreeMap<Object, Object> properties = Molecule.getProperties(molecule);
                molecule.setProperties(properties);

                Set<String> allDescriptors = new HashSet<>();
                allDescriptors.add("Molecular Formula");
                allDescriptors.add("SMILES");
                allDescriptors.add("InChI");
                allDescriptors.add("InChI Key");
                allDescriptors.add("Molecular Weight");
                allDescriptors.add("XLog P");
                allDescriptors.add("TPSA");
                allDescriptors.add("Rotatable Bond Count");
                allDescriptors.add("Hydrogen Bond Acceptor Count");
                allDescriptors.add("Hydrogen Bond Donor Count");
                allDescriptors.add("Aromatic Atom Count");
                allDescriptors.add("Aromatic Bond Count");

                SDFWriter writer = new SDFWriter(out, allDescriptors);
                writer.write(molecule);
                writer.close();

            } catch (CDKException | FileNotFoundException ex) {
                Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
            }

        }

        /**
         * Saves the input string as SDfile to a specified directory or to home
         * directory, if unspecified. Input could be either SMILES or InChI.
         *
         * @param smileOrInChI
         * @param path
         * @throws CDKException
         */
        public static void saveAsSDF(String smileOrInChI, String path) throws CDKException {

            boolean isInChI = false;

            if (smileOrInChI.startsWith("InChI")) {
                isInChI = true;
                IAtomContainer mol;
                try {
                    mol = loadMoleculeFromInChI(smileOrInChI);
                    saveAsSDF(mol, path);
                } catch (NoSuchAtomTypeException | CloneNotSupportedException | IOException ex) {
                    Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
                }

            }

            if (isInChI == false) {
                IAtomContainer mol = loadMoleculeFromSMILES(smileOrInChI);
                saveAsSDF(mol, path);
            }

        }

        /**
         * Saves the molecule object as MDL MOLfile to a specified directory or
         * to home directory, if unspecified.
         *
         * @param molecule
         * @param path
         */
        public static void saveAsMOLFile(IAtomContainer molecule, String path) {
            if (path == null) {
                path = System.getProperty("user.home") + File.separator + "output.mol";
            }

            try {
                //String outFile = path.concat(path);
                MDLV2000Writer writer = new MDLV2000Writer(new FileWriter(new File(path)));
                writer.writeMolecule(molecule);
                writer.close();
            } catch (Exception ex) {
                System.err.println("Exception" + ex.getMessage());
            }
        }

        /**
         * Saves the input string as MDL MOLfile to a specified directory or to
         * home directory, if unspecified. Input could be either SMILES or
         * InChI.
         *
         * @param smileOrInChI
         * @param path
         * @throws org.openscience.cdk.exception.CDKException
         */
        public static void saveAsMOLFile(String smileOrInChI, String path) throws CDKException {

            boolean isInChI = false;

            if (smileOrInChI.startsWith("InChI")) {
                isInChI = true;
                IAtomContainer mol;
                try {
                    mol = loadMoleculeFromInChI(smileOrInChI);
                    saveAsMOLFile(mol, path);
                } catch (NoSuchAtomTypeException | CloneNotSupportedException | IOException ex) {
                    Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
                }

            }

            if (isInChI == false) {
                IAtomContainer mol = loadMoleculeFromSMILES(smileOrInChI);
                saveAsMOLFile(mol, path);
            }

        }

        /**
         * Saves the molecule as a Chemical Markup Language (CML) file to the
         * specified path or to home directory, if unspecified.
         *
         * @param mol
         * @param path
         * @throws IOException
         * @throws CDKException
         */
        public static void saveAsCML(IAtomContainer mol, String path) throws IOException, CDKException {
            if (path == null) {
                path = System.getProperty("user.home").concat("/output.cml");
            } else {
                path = path.concat("/molecule.cml");
            }

            FileWriter output = new FileWriter(path);
            CMLWriter cmlwriter = new CMLWriter(output);
            cmlwriter.write(mol);
            cmlwriter.close();
        }

        /**
         * Saves the given molecule as PNG image file to directory specified or
         * to home directory, if unspecified. Default dimensions of the rendered
         * image is 500 x 500.
         *
         * @param molecule
         * @param path
         * @throws CDKException
         * @throws FileNotFoundException
         * @throws IOException
         */
        public static void saveAsPNG(IAtomContainer molecule, String path) throws CDKException, FileNotFoundException, IOException {

            if (path == null) {
                path = System.getProperty("user.home");
            }
            OutputStream out = new FileOutputStream(path.concat("/output.png"));

            int height = 500;
            int width = 500;

            BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

            StructureDiagramGenerator sdg = new StructureDiagramGenerator();
            sdg.setMolecule(molecule);
            sdg.generateCoordinates();
            molecule = sdg.getMolecule();

            List generators = Arrays.asList(new BasicBondGenerator(),
                    new BasicAtomGenerator(),
                    new BasicSceneGenerator()
            );

            AtomContainerRenderer renderer = new AtomContainerRenderer(generators, new AWTFontManager());
            renderer.getRenderer2DModel().set(BasicAtomGenerator.ColorByType.class, true);
            Graphics2D g2 = image.createGraphics();

            g2.fillRect(0, 0, width, height);
            renderer.paint(molecule, new AWTDrawVisitor(g2),
                    new Rectangle2D.Double(0, 0, width, height), true);
            g2.dispose();

            ImageIO.write(image, "PNG", out);

        }

        /**
         * Loads the input SMILES as an AtomContainer.
         *
         * @param smiles
         * @return
         */
        public static IAtomContainer loadMoleculeFromSMILES(String smiles) {
            try {
                SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
                //parser.setPreservingAromaticity(false);
                IAtomContainer mol = parser.parseSmiles(smiles);
                StructureDiagramGenerator gen = new StructureDiagramGenerator();
                gen.setMolecule(mol);
                gen.generateCoordinates();

                return gen.getMolecule();
            } catch (InvalidSmilesException ex) {
                System.err.println("InvalidSmilesException" + ex.getMessage());
            } catch (CDKException ex) {
                Logger.getLogger(CDKUtil.class.getName()).log(Level.SEVERE, null, ex);
            }
            return null;
        }

        /**
         * Loads the input InChI as an AtomContainer.
         *
         * @param inchi
         * @return
         * @throws CDKException
         * @throws NoSuchAtomTypeException
         * @throws CloneNotSupportedException
         * @throws IOException
         */
        public static IAtomContainer loadMoleculeFromInChI(String inchi) throws CDKException, NoSuchAtomTypeException, CloneNotSupportedException, IOException {

            InChIToStructure intostruct = InChIGeneratorFactory.getInstance().getInChIToStructure(inchi, SilentChemObjectBuilder.getInstance());
            if (intostruct.getReturnStatus() != INCHI_RET.OKAY) {
                // Structure generation failed
                throw new CDKException("Structure generation failed: " + intostruct.getReturnStatus().toString()
                        + " [" + intostruct.getMessage() + "]");
            }
            IAtomContainer mol = intostruct.getAtomContainer();
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            return gen.getMolecule();
        }

        private IO() {

        }

    }

    /**
     * Provides few example CDK AtomContainers (molecules). Methods were
     * directly taken from CDK MoleculeFactory but updated to generate 2D
     * coordinates.
     */
    public static class Factory {

        public static IAtomContainer getAlphaPinene() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9
            mol.addAtom(new Atom("C")); // 10

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6
            mol.addBond(0, 6, IBond.Order.SINGLE); // 7
            mol.addBond(3, 7, IBond.Order.SINGLE); // 8
            mol.addBond(5, 7, IBond.Order.SINGLE); // 9
            mol.addBond(7, 8, IBond.Order.SINGLE); // 10
            mol.addBond(7, 9, IBond.Order.SINGLE); // 11
            configureAtoms(mol);
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Generate an Alkane (chain of carbons with no hydrogens) of a given
         * length.
         *
         * <p>
         * This method was written by Stephen Tomkinson.
         *
         * @param chainLength The number of carbon atoms to have in the chain.
         * @return A molecule containing a bonded chain of carbons.
         *
         * @cdk.created 2003-08-15
         */
        public static IAtomContainer getAlkane(int chainLength) {
            IAtomContainer currentChain = new AtomContainer();

            //Add the initial atom
            currentChain.addAtom(new Atom("C"));

            //Add further atoms and bonds as needed, a pair at a time.
            for (int atomCount = 1; atomCount < chainLength; atomCount++) {
                currentChain.addAtom(new Atom("C"));
                currentChain.addBond(atomCount, atomCount - 1, IBond.Order.SINGLE);
            }

            return currentChain;
        }

        public static IAtomContainer getEthylCyclohexane() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6
            mol.addBond(0, 6, IBond.Order.SINGLE); // 7
            mol.addBond(6, 7, IBond.Order.SINGLE); // 8
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns cyclohexene without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C6H10/c1-2-4-6-5-3-1/h1-2H,3-6H2
         */
        public static IAtomContainer getCyclohexene() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.DOUBLE); // 6
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns cyclohexane without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2
         */
        public static IAtomContainer getCyclohexane() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns cyclopentane without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C5H10/c1-2-4-5-3-1/h1-5H2
         */
        public static IAtomContainer getCyclopentane() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.SINGLE); // 5
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns cyclobutane without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C4H8/c1-2-4-3-1/h1-4H2
         */
        public static IAtomContainer getCyclobutane() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 0, IBond.Order.SINGLE); // 4
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns cyclobutadiene without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C4H4/c1-2-4-3-1/h1-4H
         */
        public static IAtomContainer getCyclobutadiene() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.DOUBLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 0, IBond.Order.DOUBLE); // 4
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getPropylCycloPropane() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 4
            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 0, IBond.Order.SINGLE); // 3
            mol.addBond(2, 3, IBond.Order.SINGLE); // 4
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 4

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns biphenyl without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C12H10/c1-3-7-11(8-4-1)12-9-5-2-6-10-12/h1-10H
         */
        public static IAtomContainer getBiphenyl() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1		
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9
            mol.addAtom(new Atom("C")); // 10
            mol.addAtom(new Atom("C")); // 11

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6

            mol.addBond(0, 6, IBond.Order.SINGLE); // 7
            mol.addBond(6, 7, IBond.Order.SINGLE); // 8
            mol.addBond(7, 8, IBond.Order.DOUBLE); // 5
            mol.addBond(8, 9, IBond.Order.SINGLE); // 6
            mol.addBond(9, 10, IBond.Order.DOUBLE); // 7
            mol.addBond(10, 11, IBond.Order.SINGLE); // 8
            mol.addBond(11, 6, IBond.Order.DOUBLE); // 5
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getPhenylEthylBenzene() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1		
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9
            mol.addAtom(new Atom("C")); // 10
            mol.addAtom(new Atom("C")); // 11
            mol.addAtom(new Atom("C")); // 12
            mol.addAtom(new Atom("C")); // 13

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6

            mol.addBond(0, 6, IBond.Order.SINGLE); // 7
            mol.addBond(6, 7, IBond.Order.SINGLE); // 8
            mol.addBond(7, 8, IBond.Order.SINGLE); // 5
            mol.addBond(8, 9, IBond.Order.SINGLE); // 6
            mol.addBond(9, 10, IBond.Order.DOUBLE); // 7
            mol.addBond(10, 11, IBond.Order.SINGLE); // 8
            mol.addBond(11, 12, IBond.Order.DOUBLE); // 5
            mol.addBond(12, 13, IBond.Order.SINGLE);
            mol.addBond(13, 8, IBond.Order.DOUBLE); // 5
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getPhenylAmine() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1		
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("N")); // 6

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6

            mol.addBond(0, 6, IBond.Order.SINGLE); // 7
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /* build a molecule from 4 condensed triangles */
        public static IAtomContainer get4x3CondensedRings() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 0, IBond.Order.SINGLE); // 3
            mol.addBond(2, 3, IBond.Order.SINGLE); // 4
            mol.addBond(1, 3, IBond.Order.SINGLE); // 5
            mol.addBond(3, 4, IBond.Order.SINGLE); // 6
            mol.addBond(4, 2, IBond.Order.SINGLE); // 7
            mol.addBond(4, 5, IBond.Order.SINGLE); // 8
            mol.addBond(5, 3, IBond.Order.SINGLE); // 9

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getSpiroRings() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 6, IBond.Order.SINGLE); // 6
            mol.addBond(6, 0, IBond.Order.SINGLE); // 7
            mol.addBond(6, 7, IBond.Order.SINGLE); // 8
            mol.addBond(7, 8, IBond.Order.SINGLE); // 9
            mol.addBond(8, 9, IBond.Order.SINGLE); // 10
            mol.addBond(9, 6, IBond.Order.SINGLE); // 11
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getBicycloRings() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6
            mol.addBond(6, 0, IBond.Order.SINGLE); // 7
            mol.addBond(6, 7, IBond.Order.SINGLE); // 8
            mol.addBond(7, 3, IBond.Order.SINGLE); // 9
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getFusedRings() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6
            mol.addBond(5, 6, IBond.Order.SINGLE); // 7
            mol.addBond(6, 7, IBond.Order.SINGLE); // 8
            mol.addBond(7, 4, IBond.Order.SINGLE); // 9
            mol.addBond(8, 0, IBond.Order.SINGLE); // 10
            mol.addBond(9, 1, IBond.Order.SINGLE); // 11
            mol.addBond(9, 8, IBond.Order.SINGLE); // 11
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getMethylDecaline() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9
            mol.addAtom(new Atom("C")); // 10

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6
            mol.addBond(5, 6, IBond.Order.SINGLE); // 7
            mol.addBond(6, 7, IBond.Order.SINGLE); // 8RingSet
            mol.addBond(7, 8, IBond.Order.SINGLE); // 9
            mol.addBond(8, 9, IBond.Order.SINGLE); // 10
            mol.addBond(9, 0, IBond.Order.SINGLE); // 11
            mol.addBond(5, 10, IBond.Order.SINGLE); // 12
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;

        }

        public static IAtomContainer getEthylPropylPhenantren() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9
            mol.addAtom(new Atom("C")); // 10
            mol.addAtom(new Atom("C")); // 11
            mol.addAtom(new Atom("C")); // 12
            mol.addAtom(new Atom("C")); // 13
            mol.addAtom(new Atom("C")); // 14
            mol.addAtom(new Atom("C")); // 15
            mol.addAtom(new Atom("C")); // 16
            mol.addAtom(new Atom("C")); // 17
            mol.addAtom(new Atom("C")); // 18

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.DOUBLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.DOUBLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 6, IBond.Order.DOUBLE); // 6
            mol.addBond(6, 7, IBond.Order.SINGLE); // 8
            mol.addBond(7, 8, IBond.Order.DOUBLE); // 9
            mol.addBond(8, 9, IBond.Order.SINGLE); // 10
            mol.addBond(9, 0, IBond.Order.DOUBLE); // 11		
            mol.addBond(9, 4, IBond.Order.SINGLE); // 12
            mol.addBond(8, 10, IBond.Order.SINGLE); // 12
            mol.addBond(10, 11, IBond.Order.DOUBLE); // 12
            mol.addBond(11, 12, IBond.Order.SINGLE); // 12
            mol.addBond(12, 13, IBond.Order.DOUBLE); // 12
            mol.addBond(13, 7, IBond.Order.SINGLE); // 12
            mol.addBond(3, 14, IBond.Order.SINGLE); // 12
            mol.addBond(14, 15, IBond.Order.SINGLE); // 12
            mol.addBond(12, 16, IBond.Order.SINGLE); // 12		
            mol.addBond(16, 17, IBond.Order.SINGLE); // 12
            mol.addBond(17, 18, IBond.Order.SINGLE); // 12	
            configureAtoms(mol);
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getSteran() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9
            mol.addAtom(new Atom("C")); // 10
            mol.addAtom(new Atom("C")); // 11
            mol.addAtom(new Atom("C")); // 12
            mol.addAtom(new Atom("C")); // 13
            mol.addAtom(new Atom("C")); // 14
            mol.addAtom(new Atom("C")); // 15
            mol.addAtom(new Atom("C")); // 16

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 6, IBond.Order.SINGLE); // 6
            mol.addBond(6, 7, IBond.Order.SINGLE); // 8
            mol.addBond(7, 8, IBond.Order.SINGLE); // 9
            mol.addBond(8, 9, IBond.Order.SINGLE); // 10
            mol.addBond(9, 0, IBond.Order.SINGLE); // 11		
            mol.addBond(9, 4, IBond.Order.SINGLE); // 12
            mol.addBond(8, 10, IBond.Order.SINGLE); // 13
            mol.addBond(10, 11, IBond.Order.SINGLE); // 14
            mol.addBond(11, 12, IBond.Order.SINGLE); // 15
            mol.addBond(12, 13, IBond.Order.SINGLE); // 16
            mol.addBond(13, 7, IBond.Order.SINGLE); // 17
            mol.addBond(13, 14, IBond.Order.SINGLE); // 18
            mol.addBond(14, 15, IBond.Order.SINGLE); // 19
            mol.addBond(15, 16, IBond.Order.SINGLE); // 20
            mol.addBond(16, 12, IBond.Order.SINGLE); // 21

            configureAtoms(mol);
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns azulene without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C10H8/c1-2-5-9-7-4-8-10(9)6-3-1/h1-8H
         */
        public static IAtomContainer getAzulene() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 6, IBond.Order.SINGLE); // 6
            mol.addBond(6, 7, IBond.Order.DOUBLE); // 8
            mol.addBond(7, 8, IBond.Order.SINGLE); // 9
            mol.addBond(8, 9, IBond.Order.DOUBLE); // 10
            mol.addBond(9, 5, IBond.Order.SINGLE); // 11
            mol.addBond(9, 0, IBond.Order.SINGLE); // 12

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns indole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C8H7N/c1-2-4-8-7(3-1)5-6-9-8/h1-6,9H
         */
        public static IAtomContainer getIndole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("N")); // 8

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 6, IBond.Order.SINGLE); // 6
            mol.addBond(6, 7, IBond.Order.DOUBLE); // 8
            mol.addBond(7, 8, IBond.Order.SINGLE); // 9
            mol.addBond(0, 5, IBond.Order.SINGLE); // 11
            mol.addBond(8, 0, IBond.Order.SINGLE); // 12

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns pyrrole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C4H5N/c1-2-4-5-3-1/h1-5H
         */
        public static IAtomContainer getPyrrole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns pyrrole anion without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C4H4N/c1-2-4-5-3-1/h1-4H/q-1
         */
        public static IAtomContainer getPyrroleAnion() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            IAtom nitrogenAnion = new Atom("N");
            nitrogenAnion.setFormalCharge(-1);
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(nitrogenAnion); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns imidazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C3H4N2/c1-2-5-3-4-1/h1-3H,(H,4,5)/f/h4H
         */
        public static IAtomContainer getImidazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("N")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns pyrazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C3H4N2/c1-2-4-5-3-1/h1-3H,(H,4,5)/f/h4H
         */
        public static IAtomContainer getPyrazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("N")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns 1,2,4-triazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C3H4N2/c1-2-4-5-3-1/h1-3H,(H,4,5)/f/h4H
         */
        public static IAtomContainer get124Triazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("N")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("N")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns 1,2,3-triazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C2H3N3/c1-2-4-5-3-1/h1-2H,(H,3,4,5)/f/h5H
         */
        public static IAtomContainer get123Triazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("N")); // 2
            mol.addAtom(new Atom("N")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns tetrazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/CH2N4/c1-2-4-5-3-1/h1H,(H,2,3,4,5)/f/h4H
         */
        public static IAtomContainer getTetrazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("N")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("N")); // 2
            mol.addAtom(new Atom("N")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns Oxazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C3H3NO/c1-2-5-3-4-1/h1-3H
         */
        public static IAtomContainer getOxazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("O")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("N")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns Isoxazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C3H3NO/c1-2-4-5-3-1/h1-3H
         */
        public static IAtomContainer getIsoxazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("O")); // 1
            mol.addAtom(new Atom("N")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns isothiazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C3H3NS/c1-2-4-5-3-1/h1-3H
         */
        public static IAtomContainer getIsothiazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("S")); // 1
            mol.addAtom(new Atom("N")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns thiadiazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C2H2N2S/c1-3-4-2-5-1/h1-2H
         */
        public static IAtomContainer getThiadiazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("S")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("N")); // 3
            mol.addAtom(new Atom("N")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns oxadiazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C2H2N2O/c1-3-4-2-5-1/h1-2H
         */
        public static IAtomContainer getOxadiazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("O")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("N")); // 3
            mol.addAtom(new Atom("N")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns pyridine without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C3H3NO/c1-2-4-5-3-1/h1-3H
         */
        public static IAtomContainer getPyridine() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns pyridine oxide without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C5H5NO/c7-6-4-2-1-3-5-6/h1-5H
         */
        public static IAtomContainer getPyridineOxide() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.getAtom(1).setFormalCharge(1);
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("O")); // 6
            mol.getAtom(6).setFormalCharge(-1);

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6
            mol.addBond(1, 6, IBond.Order.SINGLE); // 7

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns pyrimidine without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C4H4N2/c1-2-5-4-6-3-1/h1-4H
         */
        public static IAtomContainer getPyrimidine() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("N")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns pyridazine without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C4H4N2/c1-2-4-6-5-3-1/h1-4H
         */
        public static IAtomContainer getPyridazine() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("N")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns triazine without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C4H4N2/c1-2-4-6-5-3-1/h1-4H
         */
        public static IAtomContainer getTriazine() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("N")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("N")); // 5

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.DOUBLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        /**
         * Returns thiazole without explicit hydrogens.
         *
         * @cdk.inchi InChI=1/C3H3NS/c1-2-5-3-4-1/h1-3H
         */
        public static IAtomContainer getThiazole() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("N")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("S")); // 3
            mol.addAtom(new Atom("C")); // 4

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.DOUBLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 0, IBond.Order.DOUBLE); // 5

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getSingleRing() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
//		mol.addAtom(new Atom("C")); // 6
//		mol.addAtom(new Atom("C")); // 7
//		mol.addAtom(new Atom("C")); // 8
//		mol.addAtom(new Atom("C")); // 9

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6
//		mol.addBond(5, 6, IBond.Order.SINGLE); // 7
//		mol.addBond(6, 7, IBond.Order.SINGLE); // 8
//		mol.addBond(7, 4, IBond.Order.SINGLE); // 9
//		mol.addBond(8, 0, IBond.Order.SINGLE); // 10
//		mol.addBond(9, 1, IBond.Order.SINGLE); // 11		

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getDiamantane() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9
            mol.addAtom(new Atom("C")); // 10
            mol.addAtom(new Atom("C")); // 11
            mol.addAtom(new Atom("C")); // 12
            mol.addAtom(new Atom("C")); // 13

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.SINGLE); // 6
            mol.addBond(5, 6, IBond.Order.SINGLE); // 7
            mol.addBond(6, 9, IBond.Order.SINGLE); // 8
            mol.addBond(1, 7, IBond.Order.SINGLE); // 9
            mol.addBond(7, 9, IBond.Order.SINGLE); // 10
            mol.addBond(3, 8, IBond.Order.SINGLE); // 11		
            mol.addBond(8, 9, IBond.Order.SINGLE); // 12
            mol.addBond(0, 10, IBond.Order.SINGLE); // 13
            mol.addBond(10, 13, IBond.Order.SINGLE); // 14
            mol.addBond(2, 11, IBond.Order.SINGLE); // 15
            mol.addBond(11, 13, IBond.Order.SINGLE); // 16
            mol.addBond(4, 12, IBond.Order.SINGLE); // 17		
            mol.addBond(12, 13, IBond.Order.SINGLE); // 18		

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getBranchedAliphatic() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("C")); // 7
            mol.addAtom(new Atom("C")); // 8
            mol.addAtom(new Atom("C")); // 9
            mol.addAtom(new Atom("C")); // 10
            mol.addAtom(new Atom("C")); // 11
            mol.addAtom(new Atom("C")); // 12
            mol.addAtom(new Atom("C")); // 13
            mol.addAtom(new Atom("C")); // 14
            mol.addAtom(new Atom("C")); // 15
            mol.addAtom(new Atom("C")); // 16
            mol.addAtom(new Atom("C")); // 17
            mol.addAtom(new Atom("C")); // 18

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(2, 6, IBond.Order.SINGLE); // 6
            mol.addBond(6, 7, IBond.Order.SINGLE); // 7
            mol.addBond(7, 8, IBond.Order.SINGLE); // 8
            mol.addBond(6, 9, IBond.Order.SINGLE); // 9
            mol.addBond(6, 10, IBond.Order.SINGLE); // 10
            mol.addBond(10, 11, IBond.Order.SINGLE); // 11
            mol.addBond(8, 12, IBond.Order.TRIPLE); // 12
            mol.addBond(12, 13, IBond.Order.SINGLE); // 13
            mol.addBond(11, 14, IBond.Order.SINGLE); // 14
            mol.addBond(9, 15, IBond.Order.SINGLE);
            mol.addBond(15, 16, IBond.Order.DOUBLE);
            mol.addBond(16, 17, IBond.Order.DOUBLE);
            mol.addBond(17, 18, IBond.Order.SINGLE);

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getBenzene() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("C")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5

            mol.addBond(0, 1, IBond.Order.SINGLE); // 1
            mol.addBond(1, 2, IBond.Order.DOUBLE); // 2
            mol.addBond(2, 3, IBond.Order.SINGLE); // 3
            mol.addBond(3, 4, IBond.Order.DOUBLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 0, IBond.Order.DOUBLE); // 6
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getQuinone() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("O")); // 0
            mol.addAtom(new Atom("C")); // 1
            mol.addAtom(new Atom("C")); // 2
            mol.addAtom(new Atom("C")); // 3
            mol.addAtom(new Atom("C")); // 4
            mol.addAtom(new Atom("C")); // 5
            mol.addAtom(new Atom("C")); // 6
            mol.addAtom(new Atom("O")); // 7

            mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
            mol.addBond(1, 2, IBond.Order.SINGLE); // 2
            mol.addBond(2, 3, IBond.Order.DOUBLE); // 3
            mol.addBond(3, 4, IBond.Order.SINGLE); // 4
            mol.addBond(4, 5, IBond.Order.SINGLE); // 5
            mol.addBond(5, 6, IBond.Order.DOUBLE); // 6
            mol.addBond(6, 1, IBond.Order.SINGLE); // 7
            mol.addBond(4, 7, IBond.Order.DOUBLE); // 8
            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        public static IAtomContainer getPiperidine() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("N"));
            mol.addAtom(new Atom("C"));
            mol.addAtom(new Atom("C"));
            mol.addAtom(new Atom("C"));
            mol.addAtom(new Atom("C"));
            mol.addAtom(new Atom("C"));
            mol.addAtom(new Atom("H"));

            mol.addBond(0, 1, IBond.Order.SINGLE);
            mol.addBond(1, 2, IBond.Order.SINGLE);
            mol.addBond(2, 3, IBond.Order.SINGLE);
            mol.addBond(3, 4, IBond.Order.SINGLE);
            mol.addBond(4, 5, IBond.Order.SINGLE);
            mol.addBond(5, 0, IBond.Order.SINGLE);

            mol.addBond(0, 6, IBond.Order.SINGLE);

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;

        }

        public static IAtomContainer getTetrahydropyran() throws CDKException {
            IAtomContainer mol = new AtomContainer();
            mol.addAtom(new Atom("O"));
            mol.addAtom(new Atom("C"));
            mol.addAtom(new Atom("C"));
            mol.addAtom(new Atom("C"));
            mol.addAtom(new Atom("C"));
            mol.addAtom(new Atom("C"));

            mol.addBond(0, 1, IBond.Order.SINGLE);
            mol.addBond(1, 2, IBond.Order.SINGLE);
            mol.addBond(2, 3, IBond.Order.SINGLE);
            mol.addBond(3, 4, IBond.Order.SINGLE);
            mol.addBond(4, 5, IBond.Order.SINGLE);
            mol.addBond(5, 0, IBond.Order.SINGLE);

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;

        }

        /**
         * @cdk.inchi
         * InChI=1/C5H5N5/c6-4-3-5(9-1-7-3)10-2-8-4/h1-2H,(H3,6,7,8,9,10)/f/h7H,6H2
         */
        public static IAtomContainer getAdenine() throws CDKException {
            IAtomContainer mol = new AtomContainer(); // Adenine
            IAtom a1 = mol.getBuilder().newInstance(IAtom.class, "C");
            a1.setPoint2d(new Point2d(21.0223, -17.2946));
            mol.addAtom(a1);
            IAtom a2 = mol.getBuilder().newInstance(IAtom.class, "C");
            a2.setPoint2d(new Point2d(21.0223, -18.8093));
            mol.addAtom(a2);
            IAtom a3 = mol.getBuilder().newInstance(IAtom.class, "C");
            a3.setPoint2d(new Point2d(22.1861, -16.6103));
            mol.addAtom(a3);
            IAtom a4 = mol.getBuilder().newInstance(IAtom.class, "N");
            a4.setPoint2d(new Point2d(19.8294, -16.8677));
            mol.addAtom(a4);
            IAtom a5 = mol.getBuilder().newInstance(IAtom.class, "N");
            a5.setPoint2d(new Point2d(22.2212, -19.5285));
            mol.addAtom(a5);
            IAtom a6 = mol.getBuilder().newInstance(IAtom.class, "N");
            a6.setPoint2d(new Point2d(19.8177, -19.2187));
            mol.addAtom(a6);
            IAtom a7 = mol.getBuilder().newInstance(IAtom.class, "N");
            a7.setPoint2d(new Point2d(23.4669, -17.3531));
            mol.addAtom(a7);
            IAtom a8 = mol.getBuilder().newInstance(IAtom.class, "N");
            a8.setPoint2d(new Point2d(22.1861, -15.2769));
            mol.addAtom(a8);
            IAtom a9 = mol.getBuilder().newInstance(IAtom.class, "C");
            a9.setPoint2d(new Point2d(18.9871, -18.0139));
            mol.addAtom(a9);
            IAtom a10 = mol.getBuilder().newInstance(IAtom.class, "C");
            a10.setPoint2d(new Point2d(23.4609, -18.8267));
            mol.addAtom(a10);
            IBond b1 = mol.getBuilder().newInstance(IBond.class, a1, a2, IBond.Order.DOUBLE);
            mol.addBond(b1);
            IBond b2 = mol.getBuilder().newInstance(IBond.class, a1, a3, IBond.Order.SINGLE);
            mol.addBond(b2);
            IBond b3 = mol.getBuilder().newInstance(IBond.class, a1, a4, IBond.Order.SINGLE);
            mol.addBond(b3);
            IBond b4 = mol.getBuilder().newInstance(IBond.class, a2, a5, IBond.Order.SINGLE);
            mol.addBond(b4);
            IBond b5 = mol.getBuilder().newInstance(IBond.class, a2, a6, IBond.Order.SINGLE);
            mol.addBond(b5);
            IBond b6 = mol.getBuilder().newInstance(IBond.class, a3, a7, IBond.Order.DOUBLE);
            mol.addBond(b6);
            IBond b7 = mol.getBuilder().newInstance(IBond.class, a3, a8, IBond.Order.SINGLE);
            mol.addBond(b7);
            IBond b8 = mol.getBuilder().newInstance(IBond.class, a4, a9, IBond.Order.DOUBLE);
            mol.addBond(b8);
            IBond b9 = mol.getBuilder().newInstance(IBond.class, a5, a10, IBond.Order.DOUBLE);
            mol.addBond(b9);
            IBond b10 = mol.getBuilder().newInstance(IBond.class, a6, a9, IBond.Order.SINGLE);
            mol.addBond(b10);
            IBond b11 = mol.getBuilder().newInstance(IBond.class, a7, a10, IBond.Order.SINGLE);
            mol.addBond(b11);

            StructureDiagramGenerator gen = new StructureDiagramGenerator();
            gen.setMolecule(mol);
            gen.generateCoordinates();

            mol = gen.getMolecule();

            return mol;
        }

        private static void configureAtoms(IAtomContainer mol) {
            try {
                Isotopes.getInstance().configureAtoms(mol);
            } catch (IOException ex) {
                System.err.println("IOException" + ex.getMessage());
            }
        }

        private Factory() {

        }
    }

}
