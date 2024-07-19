#Variant classification
#---
#Function to classify a variant as either a transition or transversion SNP.

#paramters
#ref: array of reference sequenced base pairs
#alt: array of alternative sequenced base pairs


#Function
#---
classify_variant = function(ref, alt){
	#define purines and pyrimidines
	purines <- c("A", "G")
	pyrimidines <- c("C", "T")

	#if the variant is a transition
	if (ref %in% purines && alt %in% purines || ref %in% pyrimidines && alt %in% pyrimidines){
		return("transition")

	#if the variant is a transversion
	} else if (ref %in% purines && alt %in% pyrimidines || ref %in% pyrimidines && alt %in% purines){
		return("transversion")
	} else {
		return("other")
	}
}
