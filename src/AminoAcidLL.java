class AminoAcidLL{
  private static String[] codons;
  char aminoAcid;
  //String[] codons;
  int[] counts;
  AminoAcidLL next;

  AminoAcidLL(){

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon 
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon){
    aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
    counts = new int[codons.length];
  }

  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops, 
   * if not passes the task to the next node. 
   * If there is no next node, add a new node to the list that would contain the codon. 
   */
  private void addCodon(String inCodon){
    codons[codons.length+1] =  inCodon;
  
  }



  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    return counts.length;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    int diff = totalCount() - inList.totalCount();
    return diff;
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts. 
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
    int diff = totalCount() - inList.totalCount();
    return diff;

    
  }

  /********************************************************************************************/
  /* Same ad above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList) {
    int diff = totalCount() - inList.totalCount();
    return diff;

  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
    char [] cc = new char[codons.length];
    int c = 0;
    for(String s: codons )
    {
      cc[c] = s.charAt(0);
    }

    return cc;
    //return new char[]{};
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
    int [] cc = new int[codons.length];
    int c = 0;
    for(String s: codons )
    {
      cc[c] = s.length();
    }

    return cc;


    //return new int[]{};
  }


  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    int check = 0;
    for(int i = 0; i<codons.length-1; i++) {
      for (int j = i+1; j<codons.length; j++) {
        if(codons[i].compareTo(codons[j])>0) {
          String temp = codons[i];
          codons[i] = codons[j];
          codons[j] = temp;
          check = 1;
        }
      }
    }
    if(check == 1)
      return false;
    else
      return true;
  }



  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
    AminoAcidLL ll = new AminoAcidLL();

    if (inSequence.length() != 3)
      ll.aminoAcid =  '0';
    inSequence = inSequence.toUpperCase();
    if (inSequence.charAt(0) == 'A') {
      if (inSequence.charAt(1) == 'A') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G')
          ll.aminoAcid =  'K';
        if (inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'N';
      }
      if (inSequence.charAt(1) == 'C') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G' || inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'T';
      }
      if (inSequence.charAt(1) == 'G') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G')
          ll.aminoAcid = 'R';
        if (inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'S';
      }
      if (inSequence.charAt(1) == 'U') {
        if (inSequence.charAt(2) == 'G')
          ll.aminoAcid = 'M';
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'I';
      }
    }
    if (inSequence.charAt(0) == 'C') {
      if (inSequence.charAt(1) == 'A') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G')
          ll.aminoAcid = 'Q';
        if (inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'H';
      }
      if (inSequence.charAt(1) == 'C') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G' || inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'P';
      }
      if (inSequence.charAt(1) == 'G') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G' || inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'R';
      }
      if (inSequence.charAt(1) == 'U') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G' || inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'L';
      }
    }
    if (inSequence.charAt(0) == 'G') {
      if (inSequence.charAt(1) == 'A') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G')
          ll.aminoAcid = 'E';
        if (inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'D';
      }
      if (inSequence.charAt(1) == 'C') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G' || inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'A';
      }
      if (inSequence.charAt(1) == 'G') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G' || inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'G';
      }
      if (inSequence.charAt(1) == 'U') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G' || inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'V';
      }
    }
    if (inSequence.charAt(0) == 'U') {
      if (inSequence.charAt(1) == 'A') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G')
          ll.aminoAcid = '*';// STOP
        if (inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'T';
      }
      if (inSequence.charAt(1) == 'C') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G' || inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'S';
      }
      if (inSequence.charAt(1) == 'G') {
        if (inSequence.charAt(2) == 'A')
          ll.aminoAcid = '*'; // STOP
        if (inSequence.charAt(2) == 'G')
          ll.aminoAcid = 'W';
        if (inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'C';
      }
      if (inSequence.charAt(1) == 'U') {
        if (inSequence.charAt(2) == 'A' || inSequence.charAt(2) == 'G')
          ll.aminoAcid = 'L';
        if (inSequence.charAt(2) == 'C' || inSequence.charAt(2) == 'U')
          ll.aminoAcid = 'F';
      }
    }
    // If we ended up not finding all 3 characters, return NULL
    return null;
  }


  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList) {

    return inList;
  }
}