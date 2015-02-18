schedule.str<-"- 19th February - Introduction and Workflows
- 26th February - Image Enhancement (A. Kaestner)
- 5th March - Basic Segmentation, Discrete Binary Structures
- 12th March - Advanced Segmentation
- 19th March - Machine Learning in Image Processing (M. Jäggi and A. Lucchi)
- 26th March - Analyzing Single Objects
- 2nd April -  Analyzing Complex Objects
- 16th April -  Spatial Distribution
- 23rd April -  Statistics and Reproducibility
- 30th April - Dynamic Experiments
- 7th May - Scaling Up / Big Data
- 21th May - Guest Lecture, Applications in Material Science
- 28th May - Project Presentations"
course.desc<-
  data.frame(Lecture=substring(gsub("\\s+"," ",strsplit(schedule.str,"\n")[[1]]),3),
             Description=c("Basic overview of the course, introduction to the basics of images and their acquisition, the importance of reproducibility and why workflows make sense for image processing",
                           "Overview of what techniques are available for assessing and improving the quality of images, specifically various filters, when to apply them, their side-effects, and how to apply them correctly",
                           "How to convert images into structures, starting with very basic techniques like threshold and exploring several automated techniques",
                           "More advanced techniques for extracting structures including basic clustering and classification approaches, and component labeling",
                           "Applying more advanced techniques from the field of Machine Learning to image processing segmentation and analysis like Support vector machines (SVM) and Markov Random Fields (MRF)",
                           "The analysis and characterization of single structures/objects after they have been segmented, including shape and orientation",
                           "What techniques are available to analyze more complicated objects with poorly defined 'shape' using Distance maps, Thickness maps, and Voronoi tesselation",
                           "Extracting meaningful information for a collection of objects like their spatial distribution, alignment, connectivity, and relative positioning",
                           "Making a statistical analysis from quantified image data, and establishing the precision of the metrics calculated, also more coverage of the steps to making an analysis reproducible",
                           "Performing tracking and registration in dynamic, changing systems covering object and image based methods",
                           "Performing large scale analyses on clusters and cloud-based machines and an introduction of how to work with 'big data' frameworks",
                           "Application of the course material to an actual scientific project from material science where we reproduce the results of a publication",
                           "The presentations of the student projects done in the class"
                           ),
             Applications=c(
               "Calculating the intensity for a folder full of images",
               "Removing detector noise from neutron images to distinguish different materials",
               "Identify cells from noise, background, and dust",
               "Identifying fat and ice crystals in ice cream images",
               "Training an algorithm to automatically identify cells",
               "Count cells and determine their average shape and volume",
               "Seperate clumps of cells, analyze vessel networks, trabecular bone, and other similar structures",
               "Quantify cells as being evenly spaced or tightly clustered or organized in sheets",
               "Determine if/how different a cancerous cell is from a healthly cell properly",
               "Turning a video of foam flow into metrics like speed, average deformation, and reorganization",
               "Performing large scale analyses using ETHs clusters and Amazons Cloud Resources, how to do anything with a terabytes of data",
               "",
               ""
               )
  )
cat(schedule.str)