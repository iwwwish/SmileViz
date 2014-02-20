/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.unibonn.vishal.main;

import de.unibonn.vishal.utils.ImageRenderer;
import javax.swing.ImageIcon;
import org.openscience.cdk.exception.CDKException;

/**
 *
 * @author vsiramshetty
 */
public class DrawMolecule {

    ImageRenderer renderer = new ImageRenderer();

    public ImageIcon getImageIcon(String smile, int w, int h) throws CDKException {
        ImageIcon icon = renderer.getImageIconFromSmile(w, h, smile);
        return icon;
    }

}
