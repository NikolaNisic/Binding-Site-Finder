package chipSeq;
import java.io.*;
import java.util.*;

public class ChipSeq {
	
	private static File binding = new File("/Users/nikolanisic/Desktop/GATA2_chr1.txt");
	private static File nonBinding = new File("/Users/nikolanisic/Desktop/not_GATA2_chr1.txt");
	private static File consensus = new File("/Users/nikolanisic/Desktop/consensus2.txt");
	private static File randomBinding = new File("/Users/nikolanisic/Desktop/random-seq-binding-sameLengths.txt");
	private static File randomNonBinding = new File("/Users/nikolanisic/Desktop/random-seq-nonBinding-sameLength.txt");
	
	
	public static void main(String[] args) throws FileNotFoundException {
		
		
		//Each ArrayList slot hold a String array. Slot [0] of each array holds the sequence description. Slot [1] of each array hold the sequence itself. 
		ArrayList<String[]> parsedBinding = FastaParser(binding);
		ArrayList<String[]> parsedNonBinding = FastaParser(nonBinding); 
		ArrayList<String[]> parsedRandomBindingSet = FastaParser(randomBinding);
		ArrayList<String[]> parsedRandomNonBindingSet = FastaParser(randomNonBinding);
		
		
		String[] parsedConsensus = consensusParser(consensus);       //feed this as arg[0] to consensusCounter
		System.out.println("\n\n\n");
		Map<String, Integer> nonBindingSetCounts = seqCounter(parsedNonBinding);  //feed this as arg[1] to consensusCounter
		Map<String, Integer> bindingSetCounts = seqCounter(parsedBinding);		//feed this as arg[1] to consensusCounter
		Map<String, Integer> randomNonBindingSetCounts = seqCounter(parsedRandomNonBindingSet);
		Map<String, Integer> randomBindingSetCounts = seqCounter(parsedRandomBindingSet);
		
		
		Map<String, Integer> nonBindingConsensusFinalCount = consensusCounter(parsedConsensus, nonBindingSetCounts);
		Map<String, Integer> bindingConsensusFinalCount = consensusCounter(parsedConsensus, bindingSetCounts);
		Map<String, Integer> randomNonBindingConsensusFinalCount = consensusCounter(parsedConsensus, randomNonBindingSetCounts);
		Map<String, Integer> randomBindingConsensusFinalCount = consensusCounter(parsedConsensus, randomBindingSetCounts);
		
		System.out.println(nonBindingConsensusFinalCount);
		System.out.println(bindingConsensusFinalCount);
		System.out.println(randomNonBindingConsensusFinalCount);
		System.out.println(randomBindingConsensusFinalCount);
		
		Map<String, Double> bindingZScores = zScore(bindingConsensusFinalCount, randomBindingConsensusFinalCount);
		Map<String, Double> nonBindingZScores = zScore(nonBindingConsensusFinalCount, randomNonBindingConsensusFinalCount);
		
		System.out.println("\n\n\n");
		System.out.println(bindingZScores);
		System.out.println(nonBindingZScores);
		
		maxScoreDifference(bindingZScores, nonBindingZScores);


	}
	
	public static double calcZScore(double exp, double random) {
		double zScore = ((exp - random)/Math.sqrt(random));
		return zScore;
	}
	
	public static Map<String, Double> zScore(Map<String, Integer> exp, Map<String, Integer> random) {
		Map<String, Double> zScoreMap = new HashMap<String, Double>();
		
		for(String conSeq : exp.keySet()) {
			zScoreMap.put(conSeq, calcZScore(exp.get(conSeq), random.get(conSeq)));
		}
		return zScoreMap;
		
	}
	
	public static double maxScoreDifference(Map<String, Double> binding, Map<String, Double> nonBinding) {
		Map<String, Double> differences = new HashMap<String, Double>();
		
		for(String seq : binding.keySet()) {
			differences.put(seq, (binding.get(seq) - nonBinding.get(seq)));
		}
		
		double maxDiff = (Collections.max(differences.values()));
		
		for(Map.Entry<String, Double> entry : differences.entrySet()) {
			if(entry.getValue() == maxDiff) {
				System.out.println("The highest difference in z-scores is: " + entry.getKey() + " = " + maxDiff);
			}
		}
		
		return maxDiff;
	}
	
	
	//This method takes in the file I created in Python of each possible consensus sequence 
	//and turns it into an array of Strings, each slot holding one of the 7^6 possible sequences.
	public static String[] consensusParser(File consensus) throws FileNotFoundException {
		Scanner sc = new Scanner(consensus);
		String[] conSeqs = new String[117649];
		sc.useDelimiter(",");
		
		for(int i = 0; i < conSeqs.length; i++) {
			conSeqs[i] = sc.next();
			System.out.print(conSeqs[i] + ",");
		}
		System.out.println("\n");
		sc.close();
		return conSeqs;
	}
	
	//This method parses the provided FASTA files and returns an ArrayList of String arrays.
	//Each String[] in the ArrayList has 2 slots. 
	//Slot [0] holds the description of the sequence
	//Slot [1] holds the sequence itself.
	public static ArrayList<String[]> FastaParser(File file) throws FileNotFoundException {
		Scanner sc = new Scanner(file);
		sc.useDelimiter(">");
		ArrayList<String[]> seqs = new ArrayList<String[]>();
		
		while(sc.hasNext()) {
			seqs.add(sc.next().split("\n", 2));
		}
		//Removing all whitespace in sequences
		for(String[] s : seqs) {
			s[1] = s[1].replaceAll("\\s+", "");
//			System.out.println(s[0]);
//			System.out.println(s[1] + "\n");
		}
		
		sc.close();
		return seqs;
	}
	
	public static boolean matches(String consensus, String key) {
		boolean result = false;
		for(int i = 0; i < key.length(); i++) {
			if((consensus.charAt(i) == 'A' && (key.charAt(i) == 'A' || key.charAt(i) == 'a')) || (consensus.charAt(i) == 'C' && (key.charAt(i) == 'C' || key.charAt(i) == 'c')) || (consensus.charAt(i) == 'G' && (key.charAt(i) == 'G' || key.charAt(i) == 'g')) || (consensus.charAt(i) == 'T' && (key.charAt(i) == 'T' || key.charAt(i) == 't')) || (consensus.charAt(i) == 'R' && (key.charAt(i) == 'A' || key.charAt(i) == 'G' || key.charAt(i) == 'a' || key.charAt(i) == 'g')) || (consensus.charAt(i) == 'Y' && (key.charAt(i) == 'C' || key.charAt(i) == 'T' || key.charAt(i) == 'c' || key.charAt(i) == 't')) || (consensus.charAt(i) == 'N' && (key.charAt(i) == 'A' || key.charAt(i) == 'C' || key.charAt(i) == 'G' || key.charAt(i) == 'T' || key.charAt(i) == 'a' || key.charAt(i) == 'c' || key.charAt(i) == 'g' || key.charAt(i) == 't'))) {
				result = true;
			}
			else {
				result = false;
			}
		}
		return result;
	}
	
	public static Map<String, Integer> consensusCounter(String[] consensus, Map<String, Integer> sequenceCounts) {
		Map<String, Integer> consensusCounts = new HashMap<String, Integer>();
		
		for(String cons : consensus) {
			for(String key : sequenceCounts.keySet()) {
				if(matches(cons, key)) {
					if(consensusCounts.containsKey(cons)) {
						consensusCounts.put(cons, consensusCounts.get(cons)+sequenceCounts.get(key));
					}
					else {
						consensusCounts.put(cons, sequenceCounts.get(key));
					}
				}
				else {
					continue;
				}
			}
		}
		return consensusCounts;
	}
	
	public static Map<String, Double> frequencyMap(Map<String, Integer> map) {
		Map<String, Double> frequencyMap = new HashMap<String, Double>();
		int totalMatches = 0;
		for(String key : map.keySet()) {
			totalMatches += map.get(key);
		}
		
		for(String key: map.keySet()) {
			frequencyMap.put(key, (double) map.get(key)/totalMatches);
		}
		return frequencyMap;
	}
	
	public static void incrementValue(Map<String, Integer> map, String key) {
		Integer count = map.get(key);
		map.put(key,  count +1);
	}
	
	public static Map<String, Integer> seqCounter(ArrayList<String[]> sequences) {
		Map<String, Integer> seqCounts = new HashMap<String, Integer>();
		
		for(String[] sequence : sequences) {
			for(int i = 0; i < sequence[1].length()-5; i++) {
				String currentSeq = (String.valueOf(sequence[1].charAt(i)) + String.valueOf(sequence[1].charAt(i+1)) + String.valueOf(sequence[1].charAt(i+2)) + String.valueOf(sequence[1].charAt(i+3)) + String.valueOf(sequence[1].charAt(i+4)) + String.valueOf(sequence[1].charAt(i+5)));
				if(seqCounts.containsKey(currentSeq)) {
					seqCounts.put(currentSeq, seqCounts.get(currentSeq)+1);
				}
				else {
					seqCounts.put(currentSeq, 1);
				}
			}
		}
		System.out.println(seqCounts);
		return seqCounts;
	}
	
}
