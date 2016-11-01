# Script to demonstrate rewriting b factor column in a pdb file

#' Overwrite b factors in a pdb file with new data
#'
#' This function reads in a pdb file, replaces the B (temperature) factors with
#' new values, then outputs the result as a new pdb file. B factors can be
#' specified per residue (the default), or by atom.
#' @param
#' in_file The pdb file to be modified
#' out_file The name of the modified pdb file to be output
#' residues A list of residues to write data to
#' values A list of the data you want to write to the pdb file. Should be same length as 'residues'!
#' pdb_chain The pdb chain to overwrite (defaults to "A")
#' @keywords pdb
#' @details
#' @examples
#' \dontrun{
#' residues=1:75
#' values=75:1
#' rewrite_pdb_bfactors("input.pdb", "output.pdb", residues, values, pdb_chain="B")
#' }
rewrite_pdb_bfactors = function(in_file, out_file, residues, values, pdb_chain="B"){
  require("Rpdb")

  # Read in the pdb file. Note, if a pdb file contains multiple models, this will not work as written
  pdb = read.pdb(in_file)

  for (i in 1:length(residues)){
    mask = (pdb[["atoms"]]$resid==residues[i] &
              pdb[["atoms"]]$chainid==pdb_chain)
    # Check to make sure at least some atoms have been selected
    # (all(!mask) only returns TRUE if all elements of mask are FALSE or NA)
    if(all(!mask))	warning("Residue ", residues[i],
                           ", chain ", pdb_chain,
                           " not found in file ", in_file)
    else pdb[["atoms"]]$temp[mask] = values[i]
  }

  write.pdb(pdb, out_file)
}


#### Main script ####
# Set the working directory to the loaction of your pdb and data file
setwd("/Users/charlie/Documents/MBiochem/R files")

# Read in the data you want to map onto the structure
# The simplest format is as a simple tab-separated table with columns for residue number and data value.
# If the import doesn't work, check that the input file is formatted correctly.
data_for_mapping = read.delim("4WYQ_chainB.txt", skip=1, header=FALSE, stringsAsFactors=FALSE)

res_ids = data_for_mapping[,1] # If you need to change the residue numbering to match the data to the pdb file, now is the time to do it.

data_values = data_for_mapping[,2]

####
#BIO3D CODE
####
# Write the data to a new pdb file
rewrite_pdb_bfactors("4WYQ_chain_B.pdb", "output.pdb", residues=res_ids, values=data_values, pdb_chain="B")
# require bio3d package
pdb <- read.pdb("xyz123") # read G-protein structure
bs <- binding.site(pdb)

print(bs$resnames) # residue names of identified binding site

