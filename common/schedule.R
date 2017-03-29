
schedule.str<-"- 23th February - Introduction and Workflows
- 2rd March - Image Enhancement (A. Kaestner)
- 9th March - Tutorial on Python and Jupyter
- 16th March - Basic Segmentation, Discrete Binary Structures
- 23th March - Advanced Segmentation
- 30th March - Analyzing Single Objects
- 6th April - Analyzing Complex Objects
- 13th April - Many Objects and Distributions
- 27th April -  Statistics and Reproducibility
- 4th May - Dynamic Experiments
- 11th May - Scaling Up / Big Data
- 18th May  - Guest Lecture - High Content Screening (M. Prummer) / Project Presentations
- 1st June - Guest Lecture - Big Aerial Images with Deep Learning and More Advanced Approaches (J. Montoya)
"
course.desc<-
  data.frame(Lecture=substring(gsub("\\s+"," ",strsplit(schedule.str,"\n")[[1]]),3),
             Description=c("Basic overview of the course, introduction to the basics of images and their acquisition, the importance of reproducibility and why workflows make sense for image processing",
                           "Overview of what techniques are available for assessing and improving the quality of images, specifically various filters, when to apply them, their side-effects, and how to apply them correctly",
                           "An introduction to the Python world of image analysis and the scikit projects",
                           "How to convert images into structures, starting with very basic techniques like threshold and exploring several automated techniques",
                           "More advanced techniques for extracting structures including basic clustering and classification approaches, and component labeling",
                           "The analysis and characterization of single structures/objects after they have been segmented, including shape and orientation",
                           "What techniques are available to analyze more complicated objects with poorly defined 'shape' using Distance maps, Thickness maps, and Voronoi tesselation",
                           "Extracting meaningful information for a collection of objects like their spatial distribution, alignment, connectivity, and relative positioning",
                           "Making a statistical analysis from quantified image data, and establishing the precision of the metrics calculated, also more coverage of the steps to making an analysis reproducible",
                           "Performing tracking and registration in dynamic, changing systems covering object and image based methods",
                           "Performing large scale analyses on clusters and cloud-based machines and an introduction of how to work with 'big data' frameworks",
                           "How Roche does Microscopy at Scale with High Content Screening and what the important image analysis aspects are",
                           "Applying more advanced techniques from the field of Machine Learning to image processing segmentation and analysis of aerial images specifically Support vector machines (SVM) and Markov Random Fields (MRF)"
                           ),
             Applications=c(
               "Calculating the intensity for a folder full of images",
               "Removing detector noise from neutron images to distinguish different materials",
               "Getting familiar with Python and learning how the basic scikit tools work",
               "Identify cells from noise, background, and dust",
               "Identifying fat and ice crystals in ice cream images",
               "Count cells and determine their average shape and volume",
               "Seperate clumps of cells, analyze vessel networks, trabecular bone, and other similar structures",
               "Quantify cells as being evenly spaced or tightly clustered or organized in sheets",
               "Determine if/how different a cancerous cell is from a healthly cell properly",
               "Turning a video of foam flow into metrics like speed, average deformation, and reorganization",
               "Performing large scale analyses using ETHs clusters and Amazons Cloud Resources, how to do anything with a terabytes of data",
               "Robust analysis of millions of images for making decisions about pharmaceuticals to pursue",
               "Identifying houses, streets, and cars in satellite images"
               )
  )
cat(schedule.str)