/*
 * Copyright (c) 2014. Vishal Siramshetty
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.unibonn.vishal.utils;

import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * This class is used to render images of chemical structures when corresponding
 * SMILES are provided as inputs.
 *
 * @author Vishal Siramshetty <srmshtty[at]gmail.com>
 */
public class ImageRenderer {

    /* default width and height for image renderer can be changed whenever needed */
    private int width = 300;
    private int height = 300;

    public Image getImageFromSmile(String smiles) throws InvalidSmilesException, CDKException {

        BufferedImage image = new BufferedImage(this.width, this.height, BufferedImage.TYPE_INT_RGB);
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);

        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        sdg.setMolecule(molecule);
        sdg.generateCoordinates();
        molecule = sdg.getMolecule();

        List generators = Arrays.asList(new SmoothGenerator(),
                new BasicSceneGenerator()
        );

        AtomContainerRenderer renderer = new AtomContainerRenderer(generators, new AWTFontManager());
        renderer.getRenderer2DModel().set(BasicSceneGenerator.FitToScreen.class, true);
        renderer.getRenderer2DModel().set(BasicSceneGenerator.UseAntiAliasing.class, true);
        renderer.getRenderer2DModel().set(BasicSceneGenerator.ZoomFactor.class, 2.0);
        

        Graphics2D g2 = image.createGraphics();

        g2.fillRect(0, 0, width, height);
        renderer.paint(molecule, new AWTDrawVisitor(g2),
                new Rectangle2D.Double(0, 0, width, height), true);
        g2.dispose();

        return image;
    }

    public ImageIcon getImageIconFromSmile(int w, int h, String smile) throws CDKException {
        setHeight(h);
        setWidth(w);
        Image im = getImageFromSmile(smile);
        return new ImageIcon(im);

    }

    /**
     * This method returns an array of images when provided a list of SMILES.
     *
     * @param smiles
     * @param height
     * @param width
     * @return Array containing images of molecules
     */
    public Image[] getImageArray(List<String> smiles, int height, int width) {
        this.height = height;
        this.width = width;

        Image[] allImages = new Image[smiles.size()];
        int i = 0;
        for (String s : smiles) {
            try {
                Image img = getImageFromSmile(s);
                allImages[i] = img;
                i++;
            } catch (CDKException ex) {
                Logger.getLogger(ImageRenderer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        return allImages;
    }

    public Icon getIconFromSmiles(String smiles) throws InvalidSmilesException, CDKException {

        return new ImageIcon(getImageFromSmile(smiles));
    }

    public int getHeight() {
        return height;
    }

    public void setHeight(int height) {
        this.height = height;
    }

    public int getWidth() {
        return width;
    }

    public void setWidth(int width) {
        this.width = width;
    }

}
