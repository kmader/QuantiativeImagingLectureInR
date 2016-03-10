schedule.str<-"- 25th February - Introduction and Workflows
- 3rd March - Image Enhancement (A. Kaestner)
- 10th March - Basic Segmentation, Discrete Binary Structures
- 17th March - Advanced Segmentation
- 24th March - Analyzing Single Objects
- 7th April - Analyzing Complex Objects
- 14th April - Spatial Distribution
- 21st April -  Statistics and Reproducibility
- 28th April - Dynamic Experiments
- 12th May - Scaling Up / Big Data
- 19th May - Guest Lecture - High Content Screening
- 26th May - Guest Lecture - Machine Learning / Deep Learning and More Advanced Approaches
- 2nd June - Project Presentations"
course.desc<-
  data.frame(Lecture=substring(gsub("\\s+"," ",strsplit(schedule.str,"\n")[[1]]),3),
             Description=c("Basic overview of the course, introduction to the basics of images and their acquisition, the importance of reproducibility and why workflows make sense for image processing",
                           "Overview of what techniques are available for assessing and improving the quality of images, specifically various filters, when to apply them, their side-effects, and how to apply them correctly",
                           "How to convert images into structures, starting with very basic techniques like threshold and exploring several automated techniques",
                           "More advanced techniques for extracting structures including basic clustering and classification approaches, and component labeling",
                           "The analysis and characterization of single structures/objects after they have been segmented, including shape and orientation",
                           "What techniques are available to analyze more complicated objects with poorly defined 'shape' using Distance maps, Thickness maps, and Voronoi tesselation",
                           "Extracting meaningful information for a collection of objects like their spatial distribution, alignment, connectivity, and relative positioning",
                           "Making a statistical analysis from quantified image data, and establishing the precision of the metrics calculated, also more coverage of the steps to making an analysis reproducible",
                           "Performing tracking and registration in dynamic, changing systems covering object and image based methods",
                           "Performing large scale analyses on clusters and cloud-based machines and an introduction of how to work with 'big data' frameworks",
                           "Application of the course material to an actual scientific project from material science where we reproduce the results of a publication",
                           "Applying more advanced techniques from the field of Machine Learning to image processing segmentation and analysis like Support vector machines (SVM) and Markov Random Fields (MRF)",
                           "The presentations of the student projects done in the class"
                           ),
             Applications=c(
               "Calculating the intensity for a folder full of images",
               "Removing detector noise from neutron images to distinguish different materials",
               "Measuring improvement in signal to noise, validating methods",
               "Identify cells from noise, background, and dust",
               "Identifying fat and ice crystals in ice cream images",
               "Count cells and determine their average shape and volume",
               "Seperate clumps of cells, analyze vessel networks, trabecular bone, and other similar structures",
               "Quantify cells as being evenly spaced or tightly clustered or organized in sheets",
               "Determine if/how different a cancerous cell is from a healthly cell properly",
               "Turning a video of foam flow into metrics like speed, average deformation, and reorganization",
               "Performing large scale analyses using ETHs clusters and Amazons Cloud Resources, how to do anything with a terabytes of data",
               "Training an algorithm to automatically identify cells",
               ""
               )
  )
cat(schedule.str)