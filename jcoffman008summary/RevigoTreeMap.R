#File mostly written by Revigo webtool- I changed the settings in the treemap function call to enhance the readabiltiy of the visualization 

# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
install.packages( "treemap", ,repos = "http://cran.us.r-project.org" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.
revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");

#below is sample data. running Revigo (provided list of GO groups and FDRs from results of GSEA or GOrilla) yields a file nearly identical to this with revigo.data specific to the data-- replace revigo.data in this file with those lines
revigo.data <- rbind(c("GO:0002088","lens development in camera-type eye",0.0143915156480057,38.6716203965613,0.699054864754172,0,"lens development in camera-type eye"),
c("GO:0048769","sarcomerogenesis",0.000495689402342331,4.0195421077239,0.725916210953988,0.37775432,"lens development in camera-type eye"),
c("GO:0007601","visual perception",0.073791629028695,28.6695862266508,0.741672870353471,0.45735326,"lens development in camera-type eye"),
c("GO:0031033","myosin filament organization",0.00108225519511409,3.70333480973847,0.851534781536639,0.48402861,"lens development in camera-type eye"),
c("GO:0001756","somitogenesis",0.0161553437713405,3.67778070526608,0.672001105019966,0.54195831,"lens development in camera-type eye"),
c("GO:0048856","anatomical structure development",1.74471103538446,3.11861534322943,0.816271967556063,0.59936794,"lens development in camera-type eye"),
c("GO:0003008","system process",0.783280132091312,15.7644715530925,0.737760567684946,0.60891247,"lens development in camera-type eye"),
c("GO:0030198","extracellular matrix organization",0.0589787773886983,6.47108329972234,0.868154933498638,0,"extracellular matrix organization"),
c("GO:0043062","extracellular structure organization",0.0606971673168184,6.03715731879876,0.90123924679827,0.32207512,"extracellular matrix organization"),
c("GO:0032501","multicellular organismal process",2.32918667117634,10.4698003017969,1,0,"multicellular organismal process"),
c("GO:0007186","G protein-coupled receptor signaling pathway",1.16410590767587,4.71896663275227,0.843657181458631,0.00874343,"G protein-coupled receptor signaling pathway"),
c("GO:0010817","regulation of hormone levels",0.195066172056766,3.70774392864352,0.88756450679122,0.2100876,"G protein-coupled receptor signaling pathway"),
c("GO:0007200","phospholipase C-activating G protein-coupled receptor signaling pathway",0.0148293746200747,3.20690839982342,0.869769926738577,0.36537001,"G protein-coupled receptor signaling pathway"),
c("GO:0070887","cellular response to chemical stimulus",1.88913427350192,3.16115090926274,0.850360784090351,0.4603421,"G protein-coupled receptor signaling pathway"),
c("GO:0007195","adenylate cyclase-inhibiting dopamine receptor signaling pathway",0.00354417922674766,3.0762380391713,0.832535153079676,0.46174904,"G protein-coupled receptor signaling pathway"),
c("GO:0042445","hormone metabolic process",0.117866678386967,3.6252516539899,0.844198167821411,0.53236442,"G protein-coupled receptor signaling pathway"),
c("GO:0034620","cellular response to unfolded protein",0.048168617672616,3.23882418684427,0.859102415735651,0.55720164,"G protein-coupled receptor signaling pathway"),
c("GO:0042026","protein refolding",0.0994022481497154,3.33441900898205,0.990941173325361,0.00914495,"protein refolding"),
c("GO:0009084","glutamine family amino acid biosynthetic process",0.621231004975565,4.50584540598156,0.945856152382999,0.01090262,"glutamine family amino acid biosynthetic process"),
c("GO:0019511","peptidyl-proline hydroxylation",0.0196416925678149,4.04000516167158,0.911290088289436,0.11135211,"glutamine family amino acid biosynthetic process"),
c("GO:0018126","protein hydroxylation",0.02466880925657,3.40230481407449,0.948346204116949,0.26015287,"glutamine family amino acid biosynthetic process"),
c("GO:0019471","4-hydroxyproline metabolic process",0.000173491290819816,4.35556141053216,0.958641938021571,0.34835457,"glutamine family amino acid biosynthetic process"),
c("GO:0007154","cell communication",7.96546432796,3.4145392704915,0.987310406351731,0.01488562,"cell communication"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- log10(as.numeric( as.character(stuff$value) ));
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  fontsize.labels=20,
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

