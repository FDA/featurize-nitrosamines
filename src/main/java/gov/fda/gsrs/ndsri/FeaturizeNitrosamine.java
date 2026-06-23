package gov.fda.gsrs.ndsri;

import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import gov.fda.gsrs.chem.util.AtomTest;
import gov.fda.gsrs.chem.util.BondTest;
import gov.fda.gsrs.chem.util.ChemUtil;
import gov.fda.gsrs.chem.util.GeomUtil;
import gov.nih.ncats.common.Tuple;
import gov.nih.ncats.common.stream.StreamUtil;
import gov.nih.ncats.common.util.CachedSupplier;
import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.AtomCoordinates;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Chemical;

/**
 * FeaturizeNitrosamine is a utility class and main executing class
 * for producing the specific features and formulae used in the
 * potency category calculations.
 *   
 * @author Tyler.Peryea
 *
 */
public class FeaturizeNitrosamine {
	
	public static class GLOBAL_SETTINGS{
		public static final String ARTIFACT_NAME = "Featureize-Nitrosamines-0.0.1-SNAPSHOT-jar-with-dependencies.jar";

		//If true, the uncharged carboxylic groups from salt splitting
		//will still be considered acidic. Otherwise only COOH will
		//count
		public static boolean CONSIDER_CHARGED_COO_AS_COOH = true;

		//If set to a positive number it will flag that rings greater than
		//this size will be considered acyclic for the "chain5" feature.
		//This is useful when dealing with large macrocycles.
		public static int MAX_RING_SIZE_BEFORE_ACYCLIC_FOR_CHAIN5 = -1;


		//Flag for whether to calculate additional features that
		//have occasionally been useful.
		public static boolean DO_EXTENDED_FEATURES_TOO = false;
	}

	public static class FeaturePairRegistry{
		public static FeatureScorePair ALPHA_HYDROGENS = new FeatureScorePair("Alpha-Hydrogens", "Alpha-hydrogen score");
		public static FeatureScorePair TERT_ALPHA_HYDROGENS = new FeatureScorePair("Tertiary alpha-carbon?", "Tertiary alpha-carbon score");
		public static FeatureScorePair COOH = new FeatureScorePair("Carboxylic acid anywhere in the molecule?", "Carboxylic acid anywhere in the molecule score");
		public static FeatureScorePair PYRROLIDINE = new FeatureScorePair("NNO in pyrrolidine ring?", "NNO in pyrrolidine ring score");

		public static FeatureScorePair S_IN6_RING = new FeatureScorePair("NNO in 6-membered ring with S?","NNO in 6-membered ring with S score");

		public static FeatureScorePair IN5_OR_6_RING = new FeatureScorePair("NNO in 5- or 6-membered ring (excluding pyrrolidine, 6-membered S-containing ring and morpholine)?", "NNO in 5- or 6-membered ring (excluding pyrrolidine, 6-membered S-containing ring and morpholine) score");

		public static FeatureScorePair MORPHOLINE = new FeatureScorePair("NNO in morpholine ring?","NNO in morpholine ring score");

		public static FeatureScorePair IN7_RING = new FeatureScorePair("NNO in a 7-membered ring?", "NNO in a 7-membered ring score");


		public static FeatureScorePair CHAIN5_BOTH = new FeatureScorePair("Chains of >=5-non-H atoms on both sides of NNO?", "Chains of >=5-non-H atoms on both sides of NNO score");

		public static FeatureScorePair EWG_ONE_SIDE = new FeatureScorePair("EWG on alpha-carbon on ONLY one side of NNO (excluding carboxylic acid, aryl and ketone)?","EWG on alpha-carbon on ONLY one side of NNO (excluding carboxylic acid, aryl and ketone) score");
		public static FeatureScorePair EWG_BOTH_SIDES = new FeatureScorePair("EWG on alpha-carbon on BOTH sides of NNO (excluding carboxylic acid, aryl and ketone)?","EWG on alpha-carbon on BOTH sides of NNO (excluding carboxylic acid, aryl and ketone) score");

		public static FeatureScorePair BETA_HYDROXYL_ONE_SIDE = new FeatureScorePair("Beta-hydroxyl on ONLY one side?","Beta-hydroxyl on ONLY one side score");
		public static FeatureScorePair BETA_HYDROXYL_BOTH_SIDES = new FeatureScorePair("Beta-hydroxyl on BOTH sides?" , "Beta-hydroxyl on BOTH sides score");

		public static FeatureScorePair ARYL_ALPHA = new FeatureScorePair("Aryl bonded to alpha-carbon?","Aryl bonded to alpha-carbon score");
		public static FeatureScorePair METHYL_BETA = new FeatureScorePair("Methyl group on beta-carbon?","Methyl group on beta-carbon score");


		//Legacy cases
		public static FeatureScorePair PIPERAZINE = new FeatureScorePair("NNO in piperazine?","NNO in piperazine score");
		public static FeatureScorePair ALLYL_GROUP = new FeatureScorePair("Double bond at beta carbon (allyl group)?","Double bond at beta carbon (allyl group) score");

	}

	

	static Map<String,SaltInfo> saltTable = SaltInfo.defaultSaltMap();

	public static Tuple<SaltInfo, Chemical> saltStrip(Chemical c, boolean makeAtomMap) throws Exception{

		c.aromatize();
		if(makeAtomMap){
			c.setAtomMapToPosition();
		}
		ChemUtil.simpleCleanup(c);

		SaltInfo sinfo = null;
		List<Chemical> components = c.connectedComponentsAsStream().collect(Collectors.toList());
		if(components.size()==1){
			sinfo=new SaltInfo(components.get(0).toInchi().getKey(), "<NO SALT>", "AM", "assumed AM");
		}else{
			boolean found=false;
			Map<String,Chemical> keep = new HashMap<>();
			for(Chemical cp:components){
				String inchiSalt = cp.toInchi().getKey();
				SaltInfo saltCMPinfo = saltTable.get(inchiSalt);
				if(saltCMPinfo!=null){
					if("USE".equalsIgnoreCase(saltCMPinfo.includeType)){
						sinfo=saltCMPinfo;
						c=cp;
						found=true;
						break;
					}else if("Exclude All".equalsIgnoreCase(saltCMPinfo.includeType)){
						sinfo=saltCMPinfo;
						sinfo.includeType="Explicitly Excluded";

						found=true;	
						break;
					}else if("Salt".equalsIgnoreCase(saltCMPinfo.includeType)){

						if(saltCMPinfo.subType.toLowerCase().contains("acid")){
							sinfo=saltCMPinfo;
						}else if(saltCMPinfo.subType.toLowerCase().contains("cation")){
							sinfo=saltCMPinfo;
						}
					}else if("Exclude".equalsIgnoreCase(saltCMPinfo.includeType)){
						//do nothing
					}else {
						keep.put(inchiSalt, cp);
					}
				}
				if(saltCMPinfo==null){
					keep.put(inchiSalt, cp);
				}				
			}

			if(!found){
				if(keep.size()==1){
					String ikey = keep.entrySet().iterator().next().getKey();
					Chemical cc = keep.entrySet().iterator().next().getValue();
					c=cc;
					if(sinfo==null){
						sinfo=new SaltInfo(ikey, "<NO ACID SALT>", "AM Stripped", "assumed AM + excluded things");
					}
				}else if(keep.size()==0){
					if(sinfo==null){
						sinfo=new SaltInfo("???", "<NOTHING>", "Nothing found", "all excluded");
					}
				}else if(keep.size()>1){
					c= keep.values().stream().reduce((ca,cb)->{
						ChemUtil.addChemical(ca, cb);
						return ca;
					}).get();
					if(sinfo==null){
						sinfo=new SaltInfo("???", "partially stripped", "Unclear which moiety to use", "");
					}
				}
			}

		}
		return Tuple.of(sinfo, c);
	}



	public static class FeatureResponse{
		private String name;
		private String type;
		private int count;
		private SaltInfo saltInfo;
		private Chemical chemical;
		private Map<String,String> featureSet = new LinkedHashMap<>();
		private Map<String,Integer> featureSetScore = new LinkedHashMap<>();

		public String getName() {
			return name;
		}

		public void setName(String name) {
			this.name = name;
		}

		public String getType() {
			return type;
		}

		public void setType(String type) {
			this.type = type;
		}

		public int getCount() {
			return count;
		}

		public void setCount(int count) {
			this.count = count;
		}

		public SaltInfo getSaltInfo() {
			return saltInfo;
		}

		public void setSaltInfo(SaltInfo saltInfo) {
			this.saltInfo = saltInfo;
		}

		public Chemical getChemical() {
			return chemical;
		}


		public Map<String, String> getFeatureSet() {
			return featureSet;
		}
		public Optional<String> getFeature(String key){
			return Optional.ofNullable(featureSet.get(key));
		}

		public void setFeatureSet(Map<String, String> featureSet) {
			this.featureSet = featureSet;
		}


		public FeatureResponse addFeature(String type, String value) {
			featureSet.put(type, value);
			return this;
		}
		public FeatureResponse addFeatureScore(String type, int value) {
			addFeature(type, value+"");
			featureSetScore.put(type, value);
			return this;
		}

		public FeatureResponse addFeatureAndScore(FeatureScorePairInstance pairInst){
			addFeature(pairInst.pair.featureName, pairInst.value);
			addFeatureScore(pairInst.pair.scoreName, pairInst.score);
			return this;
		}

		public int getSumOfScores(){
			return featureSetScore.values().stream().reduce(0, (a,b)->a+b);
		}

		/**
		 * This is where the actual algorithm flowchart is done
		 * @return
		 */
		public int getCategoryScore(){

			int alphaH   = featureSetScore.get(FeaturePairRegistry.ALPHA_HYDROGENS.scoreName);
			int tertCarb = featureSetScore.get(FeaturePairRegistry.TERT_ALPHA_HYDROGENS.scoreName);
			int potScore = getSumOfScores();


			if(alphaH==5 || alphaH==4 || tertCarb == 1){
				return 5;
			}else if(potScore<=1){
				return 1;
			}else if(potScore==2){
				return 2;
			}else if(potScore==3){
				return 3;
			}else{
				return 4;
			}
		}
	}


	public static class FeatureJob{
		private static OutputStream dummy = new OutputStream() {
			@Override
			public void write(int arg0) throws IOException {}
		};
		private String inputName;
		private Chemical c;
		private int inputForceNumber;
		private boolean useMap=false;
		private boolean addNitrosamine=false;
		private Consumer<Tuple<String,Chemical>> cons=(t)->{};
		private PrintStream outStream = new PrintStream(dummy);


		public FeatureJob(Chemical c) {
			this.c=c;
			inputName=c.getName();			
		}
		public FeatureJob(String name, Chemical c,int forceC, boolean useMap, boolean addNitrosamine, PrintStream pw, Consumer<Tuple<String,Chemical>> cons) {
			this.inputName=name;
			this.c=c;
			this.useMap=useMap;
			this.addNitrosamine=addNitrosamine;
			this.inputForceNumber=forceC;

			if(pw!=null) {
				this.outStream=pw;
			}
			if(cons!=null) {
				this.cons=cons;
			}
		}

		public static FeatureJob forOneNitrosamine(Chemical c){
			c.atoms().forEach(a->a.setAtomToAtomMap(0));
			List<Integer> sites = markAllNitrosamines(c);

			List<Chemical> chems = new ArrayList<>();
			for(int m:sites){
				chems.add(removeNitrosamine(c.copy(), m));
			}
			if(chems.size()!=1){
				System.out.println(chems.size());
				throw new IllegalArgumentException("Wrong number of nitrosamines for this method. Expected 1.");
			}
			FeatureJob fj = new FeatureJob(chems.get(0));
			fj.addNitrosamine=true;
			fj.useMap=true;

			return fj;
		}
	}



	public static List<FeatureResponse> fingerprintNitrosamine(FeatureJob fj) throws Exception{

		Chemical c = fj.c.copy();
		ChemUtil.simpleCleanup(c);


		//This shouldn't really be required
		if(!c.has2DCoordinates()) {
			try{
				c.generateCoordinates();
			}catch(Exception e){
				Map<Atom,Integer> charges = c.atoms().filter(at->at.getCharge()!=0)
						.collect(Collectors.toMap(at->at, at->at.getCharge()));
				if(charges.size()>0){
					charges.keySet().forEach(at->at.setCharge(0));
					c.generateCoordinates();
					charges.forEach((at,ch)->{
						at.setCharge(ch);
					});
				}else{
					throw e;
				}

			}
		}

		Chemical c2 = c;

		List<FeatureResponse> featureResponses = new ArrayList<>();


		Set<Integer> okayAtoms = new HashSet<>();

		if(fj.useMap){
			Atom[] arr=c.atoms().toArray(i->new Atom[i]);
			for(int i=0;i<arr.length;i++){
				if(arr[i].getAtomToAtomMap().orElse(0)>0){
					okayAtoms.add(i);
				}
			}
		}
		Tuple<SaltInfo,Chemical> retSalt= saltStrip(c,true);
		c=retSalt.v();
		SaltInfo sinfo = retSalt.k();
		c.aromatize();

		Chemical cFin = c;




		boolean got=false;
		String type = "NOT SECONDARY AMINE";

		List<AtomTest> atListDiMethyl = AtomTest.stream(c)
				.filter(at->at.isNitrogen())
				.filter(at->at.isNeutral())
				.filter(at->at.getNeighbors(nn->nn.isMethyl()).count()==2)
				.filter(at->at.getNeighbors(nn->nn.isCarbon()).count()==3)            //? N-N(C)C
				.filter(at->at.getNeighbors(nn->nn.isCarbonyl()).count()==0)          // no (C=O)-N(C)(C)

				.filter(at->at.getNeighbors(nn->nn.isDoubleBondToSNorP()).count()==0) // no C=S, C=P, C=N
				//.filter(at->at.getNeighbors(nn->nn.isTautomerOfDoubleBondedNSorPIgnoring(at)).count()==0) // no C=S, C=P, C=N

				.collect(Collectors.toList());

		if(atListDiMethyl.size()>0){

			type="D. Dimethyl-Amines";
			if(atListDiMethyl.size()>1){
				type="D. Multiple Dimethyl-Amines";
			}
			int i=1;
			for(AtomTest at: atListDiMethyl){
				int aNum=at.get().getAtomToAtomMap().orElse(0)-1;

				c2.getAtom(aNum).setAtomToAtomMap(aNum+1);
				fj.cons.accept(Tuple.of(type,c2));
				//if()
				if(!fj.useMap || okayAtoms.contains(aNum)){
					int count=(fj.inputForceNumber==0)?i:fj.inputForceNumber;

					featureResponses.add(
							calculateFeatures(c2, fj.inputName,count, c, type, sinfo,  Arrays.asList(at), fj.addNitrosamine)
							);
				}
				c2.getAtom(aNum).setAtomToAtomMap(0);
				i++;
			}
			//			printStuff(smiles, name, c, type, sinfo, atListDiMethyl);
			type = "NOT SECONDARY AMINE";
			got=true;
		}
		CachedSupplier<Map<Bond,Integer>> minRings = CachedSupplier.of(()->ChemUtil.getSmallestRingSizeForEachBond(cFin, 10)
				.entrySet()
				.stream()
				.map(Tuple::of)
				.map(Tuple.kmap(ii->cFin.getBond(ii)))
				.collect(Tuple.toMap()));

		List<AtomTest> atList = AtomTest.stream(c)
				.filter(at->at.isNitrogen())
				.filter(at->{
					if(!at.hasAromaticBond()) {
						return true;
					}else {
						BondTest bt1= at.getBonds().filter(bb->bb.isAromatic()).findFirst().orElse(null);

						Integer iRing = minRings.get().get(bt1.get());

						Set<BondTest> keepSet = new HashSet<>();

						Set<BondTest> cursorBonds = new HashSet<>();
						cursorBonds.add(bt1);
						keepSet.add(bt1);
						Set<BondTest> cursorBonds2 = new HashSet<>();

						while(!cursorBonds.isEmpty()) {
							Set<BondTest> cursorBonds2ptr=cursorBonds2;
							for(BondTest bt2:cursorBonds) {
								bt2.getNeighborBonds()
								.filter(bb->!keepSet.contains(bb))
								.filter(bb->bb.isAromatic())
								.filter(bb->iRing.equals(minRings.get().get(bb.get())))
								.forEach(bb->{
									cursorBonds2ptr.add(bb);
								}) ;													
							}
							keepSet.addAll(cursorBonds2);
							cursorBonds=cursorBonds2;
							cursorBonds2 = new HashSet<>();
						}

						boolean exoCyclicDoubleBond = keepSet.stream()
								.flatMap(bb->bb.atoms())
								.distinct()
								.flatMap(att->att.getBonds())
								.filter(bb->!keepSet.contains(bb))
								.filter(bb->!bb.isAromatic())
								.filter(bb->bb.isDoubleBond())
								.findAny().isPresent();
						return exoCyclicDoubleBond;
					}

				})
				.filter(at->at.getHCount()==1)
				.filter(at->at.getNeighbors(nn->nn.isCarbonyl()).count()==0)
				.filter(at->at.getNeighbors(nn->nn.isDoubleBondToSNorP()).count()==0)
				.filter(at->{
					if(at.hasAromaticBond()) {
						if(at.getNeighbors(nn->nn.isDoubleOrAromaticBondToSNorP()).count()==0) {
							return true;
						}
						return false;
					}
					return true;
				})
				//.filter(at->at.getNeighbors(nn->nn.isTautomerOfDoubleBondedNSorPIgnoring(at)).count()==0) // no C=S, C=P, C=N
				.filter(at->at.getNeighbors(nn->nn.isCarbon()).count()==2)
				.collect(Collectors.toList());
		if(atList.size()==1){
			type = "A. Secondary Amine";
		}else if(atList.size()>1){
			type = "A. Multiple Secondary Amine";
		}else{
			atList = AtomTest.stream(c)
					.filter(at->at.isNitrogen())
					.filter(at->at.getHCount()==1)
					.filter(at->at.getNeighbors(nn->nn.isCarbonyl()).count()==0)
					.filter(at->at.getNeighbors(nn->nn.isCarbon()).count()==2)
					.filter(at->at.getNeighbors(nn->nn.isDoubleBondToSNorP()).count()==0)
					.collect(Collectors.toList());
			if(atList.size()==1){
				type = "B. Aromatic Secondary Amine";
			}else if(atList.size()>1){
				type = "B. Multiple Aromatic Secondary Amine";
			}else{
				atList = AtomTest.stream(c)
						.filter(at->at.isNitrogen())
						.filter(at->at.getHCount()==1)
						.filter(at->at.getNeighbors(nn->nn.isCarbonyl()).count()==0)
						.filter(at->at.getNeighbors(nn->nn.isCarbon()).count()==2)
						.collect(Collectors.toList());
				if(atList.size()==1){
					type = "C. Nitrosamide-like";
				}else if(atList.size()>1){
					type = "C. Multiple Nitrosamide-like";
				}
			}
		}
		System.out.printf("computed type: %s\n", type);

		if(atList.size()>0 || !got){
			int i=1;
			for(AtomTest at: atList){
				int aNum=at.get().getAtomToAtomMap().orElse(0)-1;
				c2.getAtom(aNum).setAtomToAtomMap(aNum+1);

				if(!fj.useMap || okayAtoms.contains(aNum)){
					fj.cons.accept(Tuple.of(type,c2));
					int count=(fj.inputForceNumber==0)?i:fj.inputForceNumber;
					featureResponses.add(
							calculateFeatures(c2, fj.inputName, count, c, type, sinfo,  Arrays.asList(at), fj.addNitrosamine)
							);
				}
				c2.getAtom(aNum).setAtomToAtomMap(0);

				i++;
			}
		}

		return featureResponses;
	}

	public static void removeNitrosamine(Chemical ct){
		Set<AtomTest> toRemove = new HashSet<>();

		AtomTest.stream(ct)
		.filter(at->at.isNitrosamineAtom()) //keep nitrosamines in filter
		.forEach(att->{
			toRemove.add(att); 
			toRemove.add(att.getNeighbors(p->p.isOxygen()).findFirst().get());	
		});

		toRemove.forEach(att->{
			ct.removeAtom(att.get());	
		});
	}

	/**
	 * Removes a specific indexed nitrosamine group and also
	 * removes the marks from all other atoms. This is done
	 * to signify one group as the kind of marked secondary
	 * amine needed for the default script to work
	 * 	 * 
	 * @param ct
	 * @param m
	 * @return
	 */
	public static Chemical removeNitrosamine(Chemical ct, int m){
		Set<AtomTest> toRemove = new HashSet<>();

		AtomTest.stream(ct)
		.filter(at->at.get().getAtomToAtomMap().orElse(0)==m)
		.flatMap(at->at.getNeighbors().filter(nn->nn.isNitrosamineAtom()))
		.forEach(att->{
			toRemove.add(att); 
			toRemove.add(att.getNeighbors(p->p.isOxygen()).findFirst().get());	
		});

		toRemove.forEach(att->{
			ct.removeAtom(att.get());	
		});
		ct.atoms().filter(at->at.getAtomToAtomMap().orElse(0)!=m).forEach(at->{
			at.setAtomToAtomMap(0);
		});
		return ct;
	}

	public static List<Integer> markAllNitrosamines(Chemical ct){

		return AtomTest.stream(ct)
				.filter(at->at.isNitrogen() && at.hasNeighbor(at2->at2.isNitrosamineAtom())) //keep nitrosamines in filter
				.map(att->{
					int map = att.get().getAtomIndexInParent()+1;
					att.get().setAtomToAtomMap(map);
					return map;
				})
				.collect(Collectors.toList());		
	}

	public static void nitrosate(Chemical ct){
		Atom nCenter=null;

		double avgBond = ct.bonds().mapToDouble(b->b.getBondLength()).average().orElse(1);


		Point2D.Double nPoint=null;
		Point2D.Double oPoint=null;


		for(int i=ct.getAtomCount()-1; i>=0;i--){
			Atom caN=ct.getAtom(i);
			int amap = caN.getAtomToAtomMap().orElse(0);
			if(amap>0){
				nCenter=caN;
				Set<AtomTest> nats=AtomTest.of(nCenter)
						.getNeighbors()
						.filter(cca->cca.isMethyl())
						.collect(Collectors.toSet());
				if(nats.size()==2){
					Atom chat=nats.iterator().next().get();

					ct.removeAtom(chat);
				}
				AtomTest.of(nCenter)
				.getNeighbors()
				.filter(cca->cca.isHydrogen())

				.forEach(hN->{
					ct.removeAtom(hN.get());
				});

				double pX=caN.getAtomCoordinates().getX()-caN.getNeighbors().stream().mapToDouble(nn->nn.getAtomCoordinates().getX()).average().orElse(0);
				double pY=caN.getAtomCoordinates().getY()-caN.getNeighbors().stream().mapToDouble(nn->nn.getAtomCoordinates().getY()).average().orElse(0);
				double norm=1/Math.sqrt(pX*pX+pY*pY);
				pX=norm*avgBond*pX;
				pY=norm*avgBond*pY;

				boolean collision=true;

				int var = 0;

				while(collision && var<=3){
					collision=false;
					int nY=-((var/2)*2-1);
					boolean inv=(var%2)!=0;

					Line2D.Double realLine = new Line2D.Double(caN.getAtomCoordinates().getX(), caN.getAtomCoordinates().getY(), caN.getAtomCoordinates().getX()+nY*pX, caN.getAtomCoordinates().getY()+nY*pY);
					Line2D.Double metaphoricalLine = new Line2D.Double(0, 0, 1, 0);
					AffineTransform at= GeomUtil.getTransformFromLineToLine(metaphoricalLine, realLine, inv);

					nPoint = new Point2D.Double(1,0);
					oPoint = new Point2D.Double(1+0.5,Math.sqrt(3)/2);
					nPoint=(java.awt.geom.Point2D.Double) at.transform(nPoint, null);
					oPoint=(java.awt.geom.Point2D.Double) at.transform(oPoint, null);

					Point2D cpt=oPoint;
					Point2D ncpt=nPoint;

					if(ct.atoms()
							.map(ff->Tuple.of(ff,new Point2D.Double(ff.getAtomCoordinates().getX(),ff.getAtomCoordinates().getY())))
							.filter(t->t.v().distance(cpt)<avgBond*.5 || t.v().distance(ncpt)<avgBond*.5)
							.findAny()
							.isPresent()){
						collision=true;
						var++;
					}

				}


			}
		}


		Atom catN= ct.addAtom("N");

		catN.setAtomCoordinates(AtomCoordinates.valueOf(nPoint.getX(),nPoint.getY()));
		Atom catO= ct.addAtom("O");
		catO.setAtomCoordinates(AtomCoordinates.valueOf(oPoint.getX(),oPoint.getY()));
		Bond cbN = ct.addBond(catN, catO, Bond.BondType.DOUBLE);
		ct.addBond(catN, nCenter, Bond.BondType.SINGLE);

	}

	public static FeatureResponse calculateFeatures(Chemical cc, String name, int count,  Chemical c, String type, SaltInfo sinfo, 
			List<AtomTest> atList,
			boolean addNitrosamine) throws Exception{
		FeatureResponse resp = new FeatureResponse();
		Chemical ct=cc.copy();

		String hcounts= "N/A";
		int hcountsTOT= 0;
		String ringType= "N/A";
		String ring5_or_6= "N/A";
		String ring7= "N/A";
		String carboxylType= "N/A";
		String sulfonicAcidType= "N/A";
		String phosphateType= "N/A";
		String hasEthane= "N/A";
		String hasBenzyl= "N/A";
		int hasBetaHydroxylCount= 0 ;
		String hasAlphaSP3NoHydrogens= "N/A";
		int methylCount=0;

		String isPyrrolidine="NO";
		String isPyrrolidineSP2="NO";
		String isPiperidine="NO";
		String isPiperidineSP2="NO";
		String isPiperazine = "NO";
		String isPiperazineSP2 = "NO";
		String isMorpholine= "NO";
		String isMorpholineSP2 = "NO";
		String isAzepane = "NO";
		String isAzepaneSP2 = "NO";
		String countOfEthoxy = "N/A"; 

		boolean[] tooSmallRing=new boolean[]{false};

		String chain5="N/A";
		String ringS6="NO";


		String betaAmideEWG= "N/A";
		String betaEsterEWG= "N/A";
		String betaImineEWG= "N/A";

		Map<String,String> lookup = Arrays.stream(("0,0	5\r\n" + 
				"0,1	4\r\n" + 
				"0,2	3\r\n" + 
				"0,3	2\r\n" + 
				"1,1	4\r\n" + 
				"1,2	3\r\n" + 
				"1,3	3\r\n" + 
				"2,2	1\r\n" + 
				"3,3	1\r\n" + 
				"2,3	1\r\n").split("\n")).map(ss->ss.trim().split("\t")).collect(Collectors.toMap(k->k[0], v->v[1]));

		// N-[CH2]-C(=[O,N,S])-([H],[C],[F],[Cl],[Br],[I],[OC],[NH2],[NHC],[NC2])
		long ewg_betaCarbonylDeriv = 0;

		// N-[CH2]-[X]
		long ewg_alphaHalogen = 0;
		// N-[CH2]-[CH|C][=|#][CH|C]-c
		long ewg_aromaticConjBeta= 0;

		//TODO: Why is this here?
		long ewg_aromaticBeta = 0;



		long ewg_michaelAcceptorImine= 0;

		long ewg_michaelAcceptorCyan= 0;


		long alphaSO2EWG= 0;
		long nitrateEWG= 0;
		long quatNitrogenEWG = 0;
		long alphaNitrileEWG = 0;
		long triflouroMethylEWG = 0;

		String betaDoubleBondEWG= "N/A";
		String betaTripleBondEWG= "N/A";

		//		Beta amide EWG	Beta ester EWG	Beta imine EWG	 SO2 EWG on alpha carbon	Double bond at beta carbon (allyl group)	Triple bond at beta carbon (propargyl group)


		String[] excludeRing= new String[]{"N/A"};

		CachedSupplier<Map<Bond,Integer>> minRings = CachedSupplier.of(()->ChemUtil.getSmallestRingSizeForEachBond(c, 10)
				.entrySet()
				.stream()
				.map(Tuple::of)
				.map(Tuple.kmap(ii->c.getBond(ii)))
				.collect(Tuple.toMap()));

		Map<AtomTest,Integer> carbonylAromatic = AtomTest.stream(c)
				.filter(cat->cat.isAromatic())
				.filter(cat->cat.isDoubleBondedCarbon())
				.flatMap(at->at.getNeighbors(nat->nat.isOxygen() || nat.isSulfur()))
				.filter(cat->cat.hasDoubleBond())
				.filter(cat->cat.getHCount()==0)
				.collect(Collectors.toMap(at->at,at->at.get().getAtomicNumber()));

		if(carbonylAromatic.size()>0) {
			carbonylAromatic.keySet().forEach(at->{
				at.get().setAtomicNumber(6);
			});
			c.aromatize();
		}

		if(carbonylAromatic.size()>0) {
			carbonylAromatic.forEach((at,s)->{
				at.get().setAtomicNumber(s);
				at.get().setImplicitHCount(0);
			});
			//			c.aromatize();
		}





		if(atList.size()==1){
			hcounts=atList.stream()
					.map(at->(String)(at.getNeighbors(atn->atn.isCarbon())
							.map(nn->nn.getHCount())
							.sorted()
							.map(ii->ii+"" )
							.collect(Collectors.joining(","))))
					.collect(Collectors.joining(";"));

			String keyLookup= Arrays.stream(hcounts.split(",")).limit(2).collect(Collectors.joining(","));

			hcountsTOT = Integer.parseInt(lookup.get(keyLookup));



			if(atList.stream()
					.filter(nn->nn.isInRing())
					.filter(at->at.getNeighbors(atn->atn.isCarbon() && atn.isInRing()).count()==2)
					.filter(nn->{
						int minRingSize = nn.getBonds()

								.map(bt->minRings.get().get(bt.get()))
								.filter(Objects::nonNull)
								.min(Comparator.naturalOrder())
								.orElse(999);

						if(minRingSize<8){ //3,4,5,6,7
							if(minRingSize<4){
								excludeRing[0]="YES:" + minRingSize;	
								tooSmallRing[0]=true;
								return false; 
							}else{
								excludeRing[0]="NO:" + minRingSize;
								return true; //4,5,6,7
							}

						}else{
							excludeRing[0]="YES:" + minRingSize;
							return false;
						}
					})

					.findAny().isPresent()){
				ringType="YES";
			}else{
				ringType="NO";
			}
		}

		if(AtomTest.stream(c)
				.filter(at->at.isCarboxyl()
						|| (GLOBAL_SETTINGS.CONSIDER_CHARGED_COO_AS_COOH && at.isCarboxylCharged())
						)
				.count()>0){
			carboxylType="YES";
		}else{
			carboxylType="NO";
		}

		if(AtomTest.stream(c)
				.filter(at->at.isSulfoxyl()
						//TODO: figure out
						||at.isSulfoxylCharged()
						)
				.count()>0){
			sulfonicAcidType="YES";
		}else{
			sulfonicAcidType="NO";
		}
		
		

		if(AtomTest.stream(c)
				.filter(at->at.isPhosphate())
				.count()>0){
			phosphateType="YES";
		}else{
			phosphateType="NO";
		}

		if(atList.size()==1){
			hasEthane=atList.stream()
					.map(at->{
						//-[CH2]-[CH3]
						boolean isEthane=at.getNeighbors(atn->atn.isCarbon() 
								&& atn.getHCount()==2 
								&& (atn.getNeighbors(n2->n2.isCarbon() && n2.getHCount()==3).count()>0)
								)
								.findAny().isPresent();
						if(isEthane){
							return "YES";
						}else{
							return "NO";
						}
					})
					.collect(Collectors.joining(";"));

			hasBetaHydroxylCount=(atList.stream()
					.map(at->{
						return at.getNeighbors(atn->atn.isCarbon() 
								&& (atn.getNeighbors(n2->n2.isCarbon() && n2.hasOH() && !n2.isCarboxyl() && !n2.hasDoubleBond()).count()>0)
								)
								.count()+0;

					})
					.reduce((a,b)->a+b).orElse(0l).intValue());




			hasAlphaSP3NoHydrogens=atList.stream()
					.map(at->{
						boolean isAlphaSP3NoH=at.getNeighbors(atn->atn.isCarbon() 
								&& atn.getHCount()==0 
								&& atn.getBonds().count()==4
								)
								.findAny().isPresent();
						if(isAlphaSP3NoH){
							return "YES";
						}else{
							return "NO";
						}
					})
					.collect(Collectors.joining(";"));

			hasBenzyl= atList.stream()
					.map(at->{
						boolean isBenzyl=at.getNeighbors(atn->atn.isSp3Carbon()
								//&& atn.getNeighbors(n2->(!n2.isCarbon() && !n2.isHydrogen() && !n2.equals(at))).count()==0
								//&& (atn.getNeighbors(n2->n2.isCarbon() && n2.isAromaticCarbon()).count()>0
								&& (atn.getNeighbors(n2->n2.isAromatic()).count()>0
										)
								)
								.findAny().isPresent();
						if(isBenzyl){
							return "YES";
						}else{
							return "NO";
						}
					})
					.collect(Collectors.joining(";"));

			AtomTest oat= atList.get(0);

			// N-(C)(C)
			methylCount= (int) atList.stream()
					.flatMap(at->{
						return at.getNeighbors(atn->atn.isSp3Carbon());
					})
					.flatMap(at->{
						//filter out alpha carbons that are bonded to non-H and non-C
						if(at.getNeighbors(atn->!atn.equals(oat))
								.filter(atn->!(atn.isCarbon()||atn.isHydrogen()))
								.count()>0
								){
							return Stream.empty();
						}
						return at.getNeighbors(atn->atn.isSp3Carbon());
					})
					.filter(at->{
						//only secondary carbons
						return at.getHCount()==1;
					})
					.filter(at->{
						return !at.hasNeighbor(nn->nn.isTerminalHeteroAtom());
					})
					.flatMap(at->{
						return at.getNeighbors(atn->atn.isMethyl());
					})
					.count();

			countOfEthoxy = atList.stream()
					.filter(at-> at.getNeighbors(att->att.isSp3Carbon()).count()==2)
					.flatMap(at-> at.getNeighbors(att->att.isSp3Carbon())
							.flatMap(nn->nn.getNeighbors(nn2->nn2.isSp3Carbon()))
							.flatMap(nn->nn.getNeighbors(nn2->nn2.isOxygen() && !nn2.isOH() && !nn2.hasDoubleBond()))

							)
					.count() +"";
			boolean isRing;

			if(ringType.equals("YES") || tooSmallRing[0]){
				isRing=true;
			}else{
				if(!excludeRing[0].equals("N/A")){
					int minRingSize = Integer.parseInt(excludeRing[0].split(":")[1]);
					if(GLOBAL_SETTINGS.MAX_RING_SIZE_BEFORE_ACYCLIC_FOR_CHAIN5<=0){
						isRing=true;
					}else{
						if(minRingSize>GLOBAL_SETTINGS.MAX_RING_SIZE_BEFORE_ACYCLIC_FOR_CHAIN5){
							isRing=false;
						}else{
							isRing=true;
						}
					}
				}else{
					isRing=false;
				}
			}

			chain5 = atList.stream()
					.filter(att->!isRing)
					.filter(att->{
						List<AtomTest> atoms = att.getNeighbors().filter(n->n.isCarbon()).collect(Collectors.toList());
						if(atoms.size()==2){
							AtomTest n1=atoms.get(0);
							AtomTest n2=atoms.get(1);

							Set<AtomTest> set1= StreamUtil.with(Stream.of(n2)).and(att).stream().collect(Collectors.toSet());
							Set<AtomTest> set2= StreamUtil.with(Stream.of(n1)).and(att).stream().collect(Collectors.toSet());

							LinkedList<AtomTest> path1 = new LinkedList<>();
							Set<AtomTest> chain41 = new HashSet<AtomTest>();
							LinkedList<AtomTest> path2 = new LinkedList<>();
							Set<AtomTest> chain42 = new HashSet<AtomTest>();

							n1.getBreadthFirstPathsUntil((ll)->{
								if(ll.size()>=5){
									if(ll.stream().filter(aaa->aaa.isHydrogen()).findAny().isPresent()){
										return true;
									}else{
										path1.addAll(ll);
										return false;
									}
								}else{
									if(ll.size()==4){
										if(!ll.stream().filter(aaa->aaa.isHydrogen()).findAny().isPresent()){
											AtomTest head = ll.peek();
											chain41.add(head);
											if(chain41.size()>1){
												if(head.getNeighbors(att2n->chain41.contains(att2n)).count()>0){
													path1.addAll(ll);
												}
											}
										}

									}
								}
								return true;
							}, set1);

							n2.getBreadthFirstPathsUntil((ll)->{
								if(ll.size()>=5){
									if(ll.stream().filter(aaa->aaa.isHydrogen()).findAny().isPresent()){
										return true;
									}else{
										path2.addAll(ll);
										return false;
									}
								}else{
									if(ll.size()==4){
										if(!ll.stream().filter(aaa->aaa.isHydrogen()).findAny().isPresent()){
											AtomTest head = ll.peek();
											chain42.add(head);
											if(chain42.size()>1){
												if(head.getNeighbors(att2n->chain42.contains(att2n)).count()>0){
													path2.addAll(ll);
												}
											}
										}

									}
								}
								return true;
							}, set2);

							if(!path1.isEmpty() && !path2.isEmpty()){
								return true;
							}

						}
						return false;						
					})
					.findAny()
					.map(oo->"YES")
					.orElse("NO");


			betaAmideEWG = atList.stream()
					.filter(at-> at.getNeighbors(att->att.hasNeighbor(at2->at2.isCarbonylAmide())).count()>0)
					.findAny()
					.map(oo->"YES")
					.orElse("NO");

			betaEsterEWG = atList.stream()
					.filter(at-> at.getNeighbors(att->att.isCarbon() && att.hasNeighbor(at2->at2.isCarbonylEster())).count()>0)
					.findAny()
					.map(oo->"YES")
					.orElse("NO");

			betaImineEWG = atList.stream()
					.filter(at-> at.getNeighbors(att->att.isCarbon() && att.hasNeighbor(at2->at2.isImine())).count()>0)
					.findAny()
					.map(oo->"YES")
					.orElse("NO");
			alphaSO2EWG = atList.stream()
					.map(at-> at.getNeighbors(att->att.isCarbon() && att.hasNeighbor(at2->at2.isSulfone())).count())
					.collect(Collectors.summingInt(l->(int)(l-0)));

			nitrateEWG = atList.stream()
					.map(at-> at.getNeighbors(att->att.isCarbon() && att.hasNeighbor(at2->at2.isNitrate())).count())
					.collect(Collectors.summingInt(l->(int)(l-0)));

			quatNitrogenEWG = atList.stream()
					.map(at-> at.getNeighbors(att->att.isCarbon() && att.hasNeighbor(at2->at2.isQuatAmine())).count())
					.collect(Collectors.summingInt(l->(int)(l-0)));

			alphaNitrileEWG = atList.stream()
					.map(at-> at.getNeighbors(att->att.isCarbon() && att.hasNeighbor(at2->at2.isCyanide())).count())
					.collect(Collectors.summingInt(l->(int)(l-0)));


			triflouroMethylEWG = atList.stream()

					.map(at-> at.getNeighbors(att->att.isCarbon() && att.hasNeighbor(at2->at2.isCarbon() && at2.countNeighbor(at3->at3.isAtomSymbol("F"))==3)).count())
					.collect(Collectors.summingInt(l->(int)(l-0)));

			//			System.out.println(triflouroMethylEWG+"::");

			betaDoubleBondEWG = atList.stream()
					.filter(at-> at.getNeighbors(att->att.isCarbon() && !att.hasDoubleBond() && att.hasNeighbor(at2->at2.hasDoubleBondCarbon())).count()>0)
					.findAny()
					.map(oo->"YES")
					.orElse("NO");

			betaTripleBondEWG= atList.stream()
					.filter(at-> at.getNeighbors(att->att.isCarbon() && !att.hasDoubleBond() && att.hasNeighbor(at2->at2.hasCarbonTripleBond())).count()>0)
					.findAny()
					.map(oo->"YES")
					.orElse("NO");


			// N-[CH2]-C(=[O,N,S])-([H],[C],[F],[Cl],[Br],[I],[OC],[NH2],[NHC],[NC2])
			// String betaCarbonylDeriv = "N/A";
			ewg_betaCarbonylDeriv= atList.stream()
					.flatMap(at-> at.getNeighbors(att->att.isCarbon()))
					.filter(cat->cat.getHCount()==2 || cat.getHCount()==1)
					.flatMap(at-> at.getNeighbors(att->att.isCarbon())) //beta carbons
					.filter(cat->cat.hasBond(bb->bb.isDoubleBond() && bb.hasAtomKind(ak->ak.isNitrogen() || ak.isOxygen() || ak.isSulfur()))) // =[O,N,S]
					.filter(cat->cat.getHCount()==1 ||
					cat.hasBond(bb->bb.isSingleBond() && 
							bb.hasAtomKind(ak->ak.isHydrogen() || 
									(ak.isCarbon() && cat.countCarbon()==2 && !cat.isCarbonyl()) || //exclude ketone
									ak.isAtomSymbol("F") ||
									ak.isAtomSymbol("Cl") ||
									ak.isAtomSymbol("Br") ||
									ak.isAtomSymbol("I") ||
									(ak.isAtomSymbol("O") && ak.countCarbon()==2) ||
									(ak.isAtomSymbol("N") && ak.getHCount() == 2) ||
									(ak.isAtomSymbol("N") && ak.getHCount() == 1 && ak.countCarbon()==2) ||
									(ak.isAtomSymbol("N") && ak.getHCount() == 0 && ak.countCarbon()==3)
									))

							) // =[O,N,S]

					.count();

			// N-[CH2]-[X]
			ewg_alphaHalogen = atList.stream()
					.flatMap(at-> at.getNeighbors(att->att.isCarbon()))
					.filter(cat->cat.getHCount()==2 || cat.getHCount()==1)
					.flatMap(at-> at.getNeighbors(att->att.isHalogen()))
					.count();

			// N-[CH2]-[O,N,S]-c
			ewg_aromaticBeta = atList.stream()
					.flatMap(at-> at.getNeighbors(att->att.isCarbon()))
					.filter(cat->cat.getHCount()==2 || cat.getHCount()==1)
					.flatMap(at-> at.getNeighbors(att->att.isOxygen() || att.isNitrogen() || att.isSulfur())) //
					.filter(cat->cat.getNeighbors().filter(ccc->ccc.isAromaticCarbon()).findAny().isPresent()) //aromatic neighbor

					.count();



			// N-[CH2]-[CH|C][=|#][CH|C]-c
			ewg_aromaticConjBeta= atList.stream()
					.flatMap(at-> at.getNeighbors(att->att.isCarbon()))
					.filter(cat->cat.getHCount()==2 || cat.getHCount()==1)
					.flatMap(at-> at.getNeighbors(att->att.isCarbon())) //
					.filter(cat->cat.getHCount()<2)
					.filter(cat->cat.getNeighbors(att->!att.isHydrogen()).count()==2)
					.flatMap(at-> at.getBonds().filter(bb->bb.isDoubleBond()||bb.isTripleBond()).map(bb->bb.get().getOtherAtom(at.get())))
					.map(ca->AtomTest.of(ca))
					.filter(cat->cat.getNeighbors(att->!att.isHydrogen()).count()==2)
					.filter(cat->cat.hasNeighbor(att->att.isAromaticCarbon()))					
					.count();


			ewg_michaelAcceptorImine= atList.stream()
					.flatMap(at-> at.getNeighbors(att->att.isCarbon()))
					.filter(cat->cat.getHCount()==2 || cat.getHCount()==1)
					.flatMap(at-> at.getNeighbors(att->att.isCarbon())) //
					.filter(cat->cat.getHCount()<2)
					.filter(cat->cat.getNeighbors(att->!att.isHydrogen()).count()==2)
					.flatMap(at-> at.getBonds().filter(bb->bb.isDoubleBond()||bb.isTripleBond()).map(bb->bb.get().getOtherAtom(at.get())))
					.map(ca->AtomTest.of(ca))
					.filter(cat->cat.getNeighbors(att->!att.isHydrogen()).count()==2)
					.filter(cat->cat.hasNeighbor(cat2->cat2.hasNeighbor(att->att.isDoubleBondToSorN() || att.isDoubleBondedOxygen())))					
					.count();

			ewg_michaelAcceptorCyan= atList.stream()
					.flatMap(at-> at.getNeighbors(att->att.isCarbon()))
					.filter(cat->cat.getHCount()==2 || cat.getHCount()==1)
					.flatMap(at-> at.getNeighbors(att->att.isCarbon())) //
					.filter(cat->cat.getHCount()<2)
					.filter(cat->cat.getNeighbors(att->!att.isHydrogen()).count()==2)
					.flatMap(at-> at.getBonds().filter(bb->bb.isDoubleBond()||bb.isTripleBond()).map(bb->bb.get().getOtherAtom(at.get())))
					.map(ca->AtomTest.of(ca))
					.filter(cat->cat.getNeighbors(att->!att.isHydrogen()).count()==2)
					.filter(cat->cat.hasNeighbor(att->att.isCyanide()))					
					.count();


			if(oat.isInRing()){

				List<Bond> cbonds=oat.getNeighbors(nn->nn.isInRing())
						.map(att->att.getBondTo(oat).get())
						.map(at->at.get())
						.collect(Collectors.toList());
				Bond toRemove = cbonds.get(0);
				Atom at2 = toRemove.getOtherAtom(oat.get());

				int oHimplicit1= at2.getImplicitHCount();
				int oHimplicit2= oat.get().getImplicitHCount();


				c.removeBond(toRemove);
				//TODO: not ideal, molwitch should do this automatically 
				at2.setImplicitHCount(oHimplicit1+1);
				oat.get().setImplicitHCount(oHimplicit2+1);


				AtomTest test1=AtomTest.of(at2);
				AtomTest test2=AtomTest.of(oat.get());

				LinkedList<AtomTest> list  =test1.getFirstPathTo(test2);
				AtomTest[] ogAtom = new AtomTest[]{test2};
				String ringTxt=list.stream()
						.map(nat->{
							String btype=nat.getBondTo(ogAtom[0]).map(b->b.getType()).orElse("?");

							String ret=btype+ nat.getAtomSymbol()+(nat.isSp3Carbon()?"P":"");
							ogAtom[0]=nat;
							return ret;
						})
						.collect(Collectors.joining(""));

				c.addBond(toRemove);
				at2.setImplicitHCount(oHimplicit1);
				oat.get().setImplicitHCount(oHimplicit2);
				ringTxt=ringTxt+BondTest.of(toRemove).getType();

				if(ringTxt.equals("-N-CP-CP-CP-CP-")){
					isPyrrolidine="YES";
				}

				if(ringTxt.matches("[-]N[-]CP*[-~]CP*[-~]CP*[-~]CP*[-]")){
					if(isPyrrolidine.equals("NO")){
						isPyrrolidineSP2="YES";
					}
				}

				if(ringTxt.equals("-N-CP-CP-CP-CP-CP-")){
					isPiperidine="YES";
				}
				if(ringTxt.matches("[-]N[-]CP[-~]C[-~]CP*[-~]CP*[-~]CP[-]") || ringTxt.matches("[-]N[-]CP[-~]C[-~]CP*[-~]C[-~]CP[-]")){
					if(isPiperidine.equals("NO")){
						isPiperidineSP2="YES";
					}
				}
				if(ringTxt.matches("[-]N[-]CP[-~]CP*[-~#]CP*[-~#]CP*[-~]CP[-]")){
					if(isPiperidine.equals("NO")){
						isPiperidineSP2="YES";
					}
				}


				if(ringTxt.equals("-N-CP-CP-N-CP-CP-")){
					isPiperazine="YES";
				}
				if(ringTxt.matches("[-]N[-]CP[-~]CP*[-~]N[-~]CP*[-~]CP[-]")){
					if(isPiperazine.equals("NO")){
						isPiperidineSP2="YES";
					}
				}

				if(ringTxt.equals("-N-CP-CP-O-CP-CP-")){
					isMorpholine="YES";
				}
				if(ringTxt.matches("[-]N[-]CP[-~]CP*[-~]O[-~]CP*[-~]CP[-]")){
					if(isMorpholine.equals("NO")){
						isMorpholineSP2="YES";
					}
				}


				if(ringTxt.equals("-N-CP-CP-CP-CP-CP-CP-")){
					isAzepane="YES";
				}
				if(ringTxt.matches("[-]N[-]CP[-~]CP*[-~]CP*[-~]CP*[-~]CP*[-~]CP[-]")){
					if(isAzepane.equals("NO")){
						isAzepaneSP2="YES";
					}
				}

				if(ringTxt.contains("S") && list.size()==6){
					ringS6="YES";
				}
			}			
		}

		ring5_or_6 = "NO";
		ring7="NO";
		if(ringType.equals("YES")){
			if(excludeRing[0].equals("NO:5") || excludeRing[0].equals("NO:6")){
				if(isMorpholine.equals("NO") && isPyrrolidine.equals("NO") && ringS6.equals("NO")){
					ring5_or_6 = "YES";
				}
			}else if(excludeRing[0].equals("NO:7")){
				ring7="YES";
			}
		}


		if(addNitrosamine){
			AtomTest at1=atList.get(0);
			int amap1=at1.get().getAtomToAtomMap().orElse(0);

			ct.atoms().forEach(at->{
				if(at.getAtomToAtomMap().orElse(-1)!=amap1){
					at.setAtomToAtomMap(0);
				}
			});
			nitrosate(ct);
		}


		long ewgCount = 
				ewg_betaCarbonylDeriv +
				ewg_alphaHalogen +
				ewg_aromaticConjBeta +
				ewg_michaelAcceptorImine +
				ewg_michaelAcceptorCyan + 
				alphaSO2EWG + 
				nitrateEWG +
				quatNitrogenEWG +
				alphaNitrileEWG +
				triflouroMethylEWG;




		resp.name=name;
		resp.chemical=ct;
		resp.count=count;
		resp.type=type;
		resp.saltInfo=sinfo;


		resp.addFeatureAndScore(FeaturePairRegistry.ALPHA_HYDROGENS.getInstance(hcounts, hcountsTOT));

		resp.addFeatureAndScore(FeaturePairRegistry.TERT_ALPHA_HYDROGENS.getInstanceYesNo(hasAlphaSP3NoHydrogens, 1));
		resp.addFeatureAndScore(FeaturePairRegistry.COOH.getInstanceYesNo(carboxylType, 3));
		resp.addFeatureAndScore(FeaturePairRegistry.PYRROLIDINE.getInstanceYesNo(isPyrrolidine, 3));

		resp.addFeatureAndScore(FeaturePairRegistry.S_IN6_RING.getInstanceYesNo(ringS6, 3));


		//TODO
		resp.addFeatureAndScore(FeaturePairRegistry.IN5_OR_6_RING.getInstanceYesNo(ring5_or_6, 2));

		resp.addFeatureAndScore(FeaturePairRegistry.MORPHOLINE.getInstanceYesNo(isMorpholine, 1));
		resp.addFeatureAndScore(FeaturePairRegistry.IN7_RING.getInstanceYesNo(ring7, 1));

		resp.addFeatureAndScore(FeaturePairRegistry.CHAIN5_BOTH.getInstanceYesNo(chain5, 1));

		resp.addFeatureAndScore(FeaturePairRegistry.EWG_ONE_SIDE.getInstanceYesNo(((ewgCount==1)?"YES":"NO"), 1));
		resp.addFeatureAndScore(FeaturePairRegistry.EWG_BOTH_SIDES.getInstanceYesNo(((ewgCount>1)?"YES":"NO"), 2));

		resp.addFeatureAndScore(FeaturePairRegistry.BETA_HYDROXYL_ONE_SIDE.getInstanceYesNo(((hasBetaHydroxylCount==1)?"YES":"NO"), 1));
		resp.addFeatureAndScore(FeaturePairRegistry.BETA_HYDROXYL_BOTH_SIDES.getInstanceYesNo(((hasBetaHydroxylCount>1)?"YES":"NO"), 2));

		resp.addFeatureAndScore(FeaturePairRegistry.ARYL_ALPHA.getInstanceYesNo(hasBenzyl, -1));
		resp.addFeatureAndScore(FeaturePairRegistry.METHYL_BETA.getInstanceYesNo(((methylCount>0)?"YES":"NO"), -1));


		if(GLOBAL_SETTINGS.DO_EXTENDED_FEATURES_TOO){
			resp.addFeatureAndScore(FeaturePairRegistry.PIPERAZINE.getInstanceYesNo(isPiperazine, 0));
			resp.addFeatureAndScore(FeaturePairRegistry.ALLYL_GROUP.getInstanceYesNo(betaDoubleBondEWG, 0));
		}
		return resp;
	}



	public static void main(String[] args) throws Exception{	

		Options options = new Options();

		Option removeNitrosamines = new Option("r", "remove-nitrosamines", false, "remove any nitrosamines which are present on input smiles before evaluating");
		Option addNitrosamines = new Option("a", "add-nitrosamines", false, "add nitrosamines to smiles on final export at the predicted site");
		Option mappedSites = new Option("pn", "predict-nitrosamines", false, "predict all nitrosamine sites based on secondary amines and dimethyl amines, otherwise only use sites that have atom maps");
		Option mapFoundNitrosamine = new Option("rm", "repress-mapping-site-generation", false, "do NOT generate atom-map nitrosamine site if none is present");

		Option printOnly = new Option("p", "print-only", false, "only print out the input structures in salt-stripped form (adding/removing nitrosamines based on settings) without generating fingerprint columns");

		Option inputFile = new Option("i", "input-file", true, "input file (tab-delimited) of SMILES to process. SMILES should be second column, 1st column will be repeated.");
		Option stdInput  = new Option("is", "std-input", false, "use std input instead of a file");
		Option addHeaders = new Option("rh", "remove-headers", false, "remove headers on export TSV");
		Option noHeaders = new Option("nh", "no-headers", false, "expect no headers on import TSV");

		Option outputFile = new Option("o", "output-file", true, "output file (tab-delimited) of processed data");
		Option stdOutput  = new Option("os", "std-output", false, "use std output instead of a file");

		Option verboseOutput  = new Option("v", "verbose-output", false, "add extra verbose columns in export");

		options.addOption(removeNitrosamines);
		options.addOption(addNitrosamines);
		options.addOption(mappedSites);
		options.addOption(mapFoundNitrosamine);
		options.addOption(printOnly);
		options.addOption(inputFile);
		options.addOption(stdInput);
		options.addOption(addHeaders);
		options.addOption(noHeaders);
		options.addOption(verboseOutput);
		

		options.addOption(outputFile);
		options.addOption(stdOutput);


		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		CommandLine cmd;

		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			formatter.printHelp("java -jar " + GLOBAL_SETTINGS.ARTIFACT_NAME + " <options>", options);
			System.exit(1);
			return;
		}

		ParsedOptions popt=new ParsedOptions();
		popt.addHeaders=!cmd.hasOption("rh");
		popt.addNitrosamines=cmd.hasOption("a");
		popt.mappedSites=!cmd.hasOption("pn");
		popt.mapFoundNitrosamine = !cmd.hasOption("rm");

		popt.printOnly=cmd.hasOption("p");

		popt.inputFile = cmd.getOptionValue("i");
		popt.stdInput = cmd.hasOption("is");

		popt.noHeaders = cmd.hasOption("nh");
		popt.outputFile = cmd.getOptionValue("o");
		popt.stdOutput = cmd.hasOption("os");
		
		popt.verboseOutput = cmd.hasOption("v");
		 
		Stream<String> inputStreamStrings = null;

		if(popt.stdInput){
			inputStreamStrings = new BufferedReader(
					new InputStreamReader(System.in, StandardCharsets.UTF_8))
					.lines();
		}else{
			try(BufferedReader red = new BufferedReader(
					new InputStreamReader(new FileInputStream(popt.inputFile), StandardCharsets.UTF_8))){
				inputStreamStrings = red
						.lines();
				run(inputStreamStrings, popt);
			}catch(Exception e){
				System.err.println("Unable to read input file, please specify a valid input file or use standard input");
				formatter.printHelp("java -jar " + GLOBAL_SETTINGS.ARTIFACT_NAME + " <options>", options);
				System.exit(1);
			}
		}

		


	}

	public static class ParsedOptions{

		boolean mapFoundNitrosamine;
		boolean removeNitrosamines;
		boolean addNitrosamines;
		boolean mappedSites;
		boolean printOnly;

		String inputFile;
		String outputFile;

		boolean addHeaders;
		boolean noHeaders;
		boolean stdInput;
		boolean stdOutput;

		boolean verboseOutput;

	}

	public static void printOnly(Stream<String> inputStream, PrintStream outPw, ParsedOptions parsedOptions) throws Exception{

		try{
			if(parsedOptions.addHeaders){
				outPw.println("Column_1" +"\t" + "Original_SMILES" + "\t" + "Modified SMILES" + "\t" + "Modified InChIKey");
			}
			inputStream.forEach(ss->{
				String[] cols = ss.split("\t");
				
				Chemical c;
				try {
					c = Chemical.parse(cols[1].trim());
				} catch (IOException e1) {
					throw new RuntimeException(e1);
				}

				if(parsedOptions.removeNitrosamines){
					removeNitrosamine(c);
				}

				Chemical cStripped;
				try {
					cStripped = saltStrip(c,false).v();
				} catch (Exception e1) {
					throw new RuntimeException(e1);
				}

				if(parsedOptions.addNitrosamines){
					nitrosate(cStripped);
				}

				try {
					outPw.println(cols[0].trim() + "\t" +cols[1].trim() + "\t"+ cStripped.toSmiles() + "\t" + cStripped.toInchi().getKey());
				} catch (Exception e) {
					e.printStackTrace();
				}
			});
		}catch(Exception e){

		}finally{
			try{
				outPw.close();
			}catch(Exception e){
				e.printStackTrace();
			}
		}
	}

	public static void setupSTDFilter(){

		PrintStream real  = System.out;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();

		System.setOut(new PrintStream(baos){
			@Override
			public void println(String s){
				if(!s.equals("Thread timeout!")){
					real.println(s);
				}
				try{
					baos.reset();
					baos.flush();
				}catch(Exception e){

				}
			}
			@Override
			public void print(String s){
				real.print(s);
				try{
					baos.reset();
					baos.flush();
				}catch(Exception e){

				}
			}
		});
	}



	public static void run(Stream<String> inputStream, ParsedOptions parsedOptions) throws Exception{
		setupSTDFilter();

		if(!parsedOptions.noHeaders){
			inputStream=inputStream.skip(1);
		}
		PrintStream outPwF=null; 
		if(parsedOptions.outputFile!=null){
			outPwF = new PrintStream(parsedOptions.outputFile);
		}else{
			if(parsedOptions.stdOutput){
				outPwF = System.out;
				outPwF.flush();
			}else{
				//ERROR
				System.err.println("NO OUTPUT FILE OR STD OUTPUT SPECIFIED, EXITING");
				System.exit(0);
			}
		}	
		PrintStream outPw=outPwF;

		if(parsedOptions.printOnly){
			printOnly(inputStream,outPw, parsedOptions);
		}else{
			if(parsedOptions.addHeaders){
				String headerLine = "Structure_Name"+ "\t" +
						"Nitrosamine Structure	NNO Instance	"
						+ (parsedOptions.verboseOutput? "Amine Category\t" : "")
						+ "Potency Category	Potency Score	"
						+ "Alpha-Hydrogens	"
						+ "Alpha-hydrogen score	"
						+ "Tertiary alpha-carbon?	"
						+ "Tertiary alpha-carbon score	"
						+ "Carboxylic acid group anywhere on molecule?	"
						+ "Carboxylic acid group anywhere on molecule score	NNO in pyrrolidine ring?	"
						+ "NNO in pyrrolidine ring score	"
						+ "NNO in 6-membered ring with S?	"
						+ "NNO in 6-membered ring with S score	"
						+ "NNO in 5- or 6-membered ring (excluding pyrrolidine, 6-membered S-containing ring and morpholine)?	"
						+ "NNO in 5- or 6-membered ring (excluding pyrrolidine, 6-membered S-containing ring and morpholine) score	"
						+ "NNO in morpholine ring?	"
						+ "NNO in morpholine ring score	"
						+ "NNO in a 7-membered ring?	"
						+ "NNO in a 7-membered ring score	"
						+ "Chains of >=5-non-H atoms on both sides of NNO?	"
						+ "Chains of >=5-non-H atoms on both sides of NNO score	"
						+ "EWG on alpha-carbon on ONLY one side of NNO (excluding carboxylic acid, aryl and ketone)?	"
						+ "EWG on alpha-carbon on ONLY one side of NNO (excluding carboxylic acid, aryl and ketone) score	"
						+ "EWG on alpha-carbon on BOTH sides of NNO (excluding carboxylic acid, aryl and ketone)?	"
						+ "EWG on alpha-carbon on BOTH sides of NNO (excluding carboxylic acid, aryl and ketone) score	"
						+ "Beta-hydroxyl on ONLY one side?	"
						+ "Beta-hydroxyl on ONLY one side score	Beta-hydroxyl on BOTH sides?	"
						+ "Beta-hydroxyl on BOTH sides score	"
						+ "Aryl bonded to alpha-carbon?	"
						+ "Aryl bonded to alpha-carbon score	"
						+ "Methyl group on beta-carbon?	"
						+ "Methyl group on beta-carbon score";

				outPw.println(headerLine);
			}

			inputStream.forEach(ss->{
				String[] cols = ss.trim().split("\t");
				if(cols.length>=2){
					Chemical c;
					List<Chemical> readList= new ArrayList<>();
					try {
						String inp=cols[1].trim();
						c = Chemical.parse(inp);
					} catch (Exception e1) {
						return;
						//throw new RuntimeException(e1);
					}

					if(parsedOptions.removeNitrosamines){
						removeNitrosamine(c);
					}

					if(parsedOptions.mappedSites && parsedOptions.mapFoundNitrosamine){
						boolean anyMapped = c.atoms()
								.filter(att->att.getAtomToAtomMap().isPresent())
								.findAny().isPresent();

						if(!anyMapped){
							List<Integer> sites = markAllNitrosamines(c);

							for(int m:sites){
								readList.add(removeNitrosamine(c.copy(), m));
							}
						}
					}
					if(readList.size()==0){
						readList.add(c);
					}

					try {
						for(int ci=0;ci<readList.size();ci++){
							Chemical c1 = readList.get(ci);

							int fnum=0;
							if(cols.length>=3){
								try{
									int input = Integer.parseInt(cols[2]);
									fnum=input;
								}catch(Exception e){

								}

							}
							if(fnum>0 && readList.size()>1 && fnum-1!=ci){
								continue;
							}
							FeatureJob fjob = new FeatureJob(cols[0],c1,fnum, parsedOptions.mappedSites, parsedOptions.addNitrosamines, outPw, (ccc)->{});
							List<FeatureResponse> resp = fingerprintNitrosamine(fjob);

							resp.forEach(frr->{
								FeatureResponse fr=frr;

								String smiles=null;
								try {
									smiles = fr.getChemical().toSmiles();
								} catch (IOException e) {
									e.printStackTrace();
								}

								fjob.outStream.print(fr.name + "\t");
								fjob.outStream.print(smiles + "\t");
								fjob.outStream.print(fr.count + "\t");
								//TODO: also export other stuff like salt info
								// if desired
								if(parsedOptions.verboseOutput){
									fjob.outStream.print(fr.getType() + "\t");
								}
								fjob.outStream.print(fr.getCategoryScore() + "\t");
								fjob.outStream.print(fr.getSumOfScores() + "\t");
								fjob.outStream.print(fr.getFeatureSet().values().stream().collect(Collectors.joining("\t")));
								fjob.outStream.println("");
								fjob.outStream.flush();


							});
						}



					} catch (Exception e) {
						e.printStackTrace();
					}					
				}else{
					//throw new RuntimeException("Expected 2 columns in input, found:" + cols.length);
				}
			});

		}
	}

	public static Optional<FeatureResponse> forMostPotentNitrosamine(Chemical c){
		c.atoms().forEach(a->a.setAtomToAtomMap(0));
		List<Integer> sites = markAllNitrosamines(c);
		if(sites.isEmpty()){
			throw new IllegalArgumentException("Wrong number of nitrosamines for this method. Expected 1 or more.");
		}

		List<Chemical> chems = new ArrayList<>();
		for(int m:sites){
			chems.add(removeNitrosamine(c.copy(), m));
		}

		return chems.stream()
				.map(cc->new FeatureJob(cc))
				.peek(fj->fj.addNitrosamine=false)
				.peek(fj->fj.useMap=true)
				.map(fj->{

					try{
						return FeaturizeNitrosamine.fingerprintNitrosamine(fj);
					}catch(Exception e){
						return new ArrayList<FeatureResponse>();
					}
				})
				.flatMap(r->r.stream())
				.filter(r->!(r.getType().startsWith("B") || r.getType().startsWith("C"))) //eliminate amides and amide-like items
				.sorted(Comparator.comparing(fr->fr.getCategoryScore()))
				.findFirst()
				.map(fr->{
					if(fr.getCategoryScore()==1){
						fr.addFeature("AI Limit (US)", "26.5 ng/day");
					}else if(fr.getCategoryScore()==2){
						fr.addFeature("AI Limit (US)", "100 ng/day");
					}else if(fr.getCategoryScore()==3){
						fr.addFeature("AI Limit (US)", "400 ng/day");
					}else if(fr.getCategoryScore()==4){
						fr.addFeature("AI Limit (US)", "1500 ng/day");
					}else if(fr.getCategoryScore()==5){
						fr.addFeature("AI Limit (US)", "1500 ng/day");
					}
					if(chems.size()==1){
						fr.setType("Single N-nitroso group");
					}else{
						fr.setType("Most potent of multiple N-nitroso groups");
					}
					return fr;
				});
	}

}