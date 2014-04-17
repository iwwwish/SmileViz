SMILES Visualizer
========

SMILES Visualizer (SmileViz) is a simple, always-on-top, Java-based desktop application that dynamically visualizes the SMILES input provided. It is very useful in interpreting SMILES step by step (atom-by-atom) and checking how the connectivity proceeds. It uses [CDK](http://sourceforge.net/apps/mediawiki/cdk/index.php?title=Main_Page)'s [AtomContainerRenderer](http://pele.farmbio.uu.se/nightly-1.3.1/cdk-javadoc-1.4.0/org/openscience/cdk/renderer/AtomContainerRenderer.html) to draw the image from SMILES. Thanks to [John May](https://plus.google.com/+JohnMay/about) for the SmoothRenderer. The light-weight visualizer differs from the commercially available structure visualizers in the following ways:

1. You need not hit enter each time you update the smile or enter a new atom (visualization is dynamically updated),
2. Any error in the smile will be displayed (dynamically),
3. You can instantly save the molecule as PNG, SDF or MOL by a simple right-click,
4. Lastly, you will not need Internet.


A demo for SMILES Visualizer can be found in my [blogpost](http://vishalkpp.blogspot.de/2014/02/simple-smiles-visualizer.html).
