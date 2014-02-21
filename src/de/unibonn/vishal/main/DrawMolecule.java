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
package de.unibonn.vishal.main;

import de.unibonn.vishal.utils.ImageRenderer;
import javax.swing.ImageIcon;
import org.openscience.cdk.exception.CDKException;

/**
 *
 * @author Vishal Siramshetty <srmshtty[at]gmail.com>
 */
public class DrawMolecule {

    ImageRenderer renderer = new ImageRenderer();

    public ImageIcon getImageIcon(String smile, int w, int h) throws CDKException {
        ImageIcon icon = renderer.getImageIconFromSmile(w, h, smile);
        return icon;
    }

}
