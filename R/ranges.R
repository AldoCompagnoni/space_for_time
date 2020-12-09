library(BIEN)

psme <- BIEN_occurrence_species(
          species = "Pseudotsuga menziesii",
          cultivated = F,
          only.new.world = F,
          all.taxonomy = T,
          native.status = T,
          observation.type = T,
          political.boundaries = T)

map('world',fill=T , col= "grey", bg="light blue")

points( cbind(Xanthium_strumarium_full$longitude,
              Xanthium_strumarium_full$latitude),
        col="red", pch=20, cex=1 ) 
points(cbind(Xanthium_strumarium_full$longitude,Xanthium_strumarium_full$latitude),col="red",pch=20,cex=1) 