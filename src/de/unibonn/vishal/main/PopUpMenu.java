/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.unibonn.vishal.main;

import de.unibonn.vishal.utils.Utility;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import org.openscience.cdk.interfaces.IAtomContainer;
import de.unibonn.vishal.utils.CDKUtil;

/**
 *
 * @author vsiramshetty
 */
public class PopUpMenu extends JPopupMenu {

    String smile;
    ImageIcon icon;

    public PopUpMenu(String smile, ImageIcon icon) {
        this.smile = smile;
        //this.label = label;
        this.icon = icon;
        initComponents();
    }

    private void initComponents() {
        JMenuItem item = new JMenuItem("save as Image");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                File fileToSave = Utility.UI.saveFile(null);
                Image img = icon.getImage();
                BufferedImage bi = new BufferedImage(img.getWidth(null), img.getHeight(null), BufferedImage.TYPE_INT_ARGB);
                Graphics2D g2 = bi.createGraphics();
                g2.drawImage(img, 0, 0, null);
                g2.dispose();
                try {
                    ImageIO.write(bi, "png", fileToSave);
                } catch (IOException ex) {
                    //Logger.getLogger(Test.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        });
        add(item);

        item = new JMenuItem("save as SDF");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                IAtomContainer molecule = CDKUtil.IO.loadMoleculeFromSMILES(smile);
                File fileToSave = Utility.UI.saveFile(null);
                CDKUtil.IO.saveAsSDFWithProperties(molecule, fileToSave.getAbsolutePath());
                //Utility.UI.showInfoMessage(null, "SD file saved to directory: \"" + System.getProperty("user.home") + "/\"");
            }
        });
        add(item);

        item = new JMenuItem("save as MOL");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                IAtomContainer molecule = CDKUtil.IO.loadMoleculeFromSMILES(smile);
                File fileToSave = Utility.UI.saveFile(null);
                CDKUtil.IO.saveAsMOLFile(molecule, fileToSave.getAbsolutePath());
                //Utility.UI.showInfoMessage(null, "MOL file saved to directory: \"" + System.getProperty("user.home") + "/\"");
            }
        });
        add(item);
    }

}
