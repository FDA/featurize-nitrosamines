package gov.nih.ncats.molwitch.renderer;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.apache.commons.io.IOUtils;
import org.junit.Ignore;
import org.junit.Test;

import gov.fda.gsrs.chem.util.AtomTest;
import gov.fda.gsrs.chem.util.BondTest;
import gov.fda.gsrs.ndsri.FeaturizeNitrosamine;
import gov.fda.gsrs.ndsri.FeaturizeNitrosamine.FeatureJob;
import gov.fda.gsrs.ndsri.FeaturizeNitrosamine.FeatureResponse;
import gov.nih.ncats.molwitch.Chemical;

import static org.junit.Assert.*;

public class FeatureParsingTests {


	private static class TestMol {
		private String molFileName;
		private boolean expectData;
		private String description;

		public TestMol(String molFileName, boolean expectData, String description) {
			this.molFileName = molFileName;
			this.expectData = expectData;
			this.description = description;
		}

		public String getMolFileName() {
			return molFileName;
		}

		public void setMolFileName(String molFileName) {
			this.molFileName = molFileName;
		}

		public boolean isExpectData() {
			return expectData;
		}

		public void setExpectData(boolean expectData) {
			this.expectData = expectData;
		}

		public String getDescription() {
			return description;
		}

		public void setDescription(String description) {
			this.description = description;
		}

	}

    @Test
    @Ignore
    public void testExportAtomMaps() throws IOException {
    	Chemical c1= Chemical.parse("[NH:1]1C2=CC=CC=C2N=C1C3=CSC=N3");
    	    	
    	String smi1= c1.toSmiles();
    	assertTrue(smi1.contains(":1"));
    }
    
    @Test
    public void testAbacavirIsSecondaryAmine() throws Exception {
    	Chemical c1= Chemical.parse("NC1=NC2=C(N=CN2[C@@H]3C[C@H](CO)C=C3)C(NC4CC4)=N1");
    	    
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	FeatureResponse resp1= resp.get(0);
    	
    	assertEquals(1,resp.size());
    	assertEquals("A. Secondary Amine" ,resp1.getType());
    	assertEquals("0,1" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.ALPHA_HYDROGENS.getFeatureName()).orElse(null));
    }
    

    @Test
    public void testIVACAFTORType() throws Exception {
    	Chemical c1= Chemical.parse("c1ccc2c(c1)c(=O)cc[nH]2");
    	    
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	FeatureResponse resp1= resp.get(0);
    	resp1.getChemical();

    	assertEquals(1,resp.size());
    	assertEquals("A. Secondary Amine" ,resp1.getType());
    	assertEquals("0,1" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.ALPHA_HYDROGENS.getFeatureName()).orElse(null));
    }
    
    @Test
    public void testFostamatibib() throws Exception {
    	Chemical c1= Chemical.parse("O.O.O.O.O.O.[Na+].[Na+].COC1=CC([NH:1]C2=NC=C(F)C(NC3=NC4=C(OC(C)(C)C(=O)N4COP([O-])([O-])=O)C=C3)=N2)=CC(OC)=C1OC");
    	    
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	assertEquals(2,resp.size());
    	FeatureResponse resp1= resp.get(0);
    	assertEquals("A. Multiple Secondary Amine" ,resp1.getType());
    	assertEquals("0,0" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.ALPHA_HYDROGENS.getFeatureName()).orElse(null));
    }


	@Test
    public void testSP3CarbonOnSmallRing() throws Exception {
    	Chemical c1= Chemical.parse("C1CCCCC1");
    	
    	assertTrue(AtomTest.stream(c1)
    	.allMatch(c->c.isSp3Carbon()));
    }
    @Test
    public void testSP3CarbonOnSmallHeteroRing() throws Exception {
    	Chemical c1= Chemical.parse("N1CCCCC1");
    	
    	assertTrue(AtomTest.stream(c1)
    			.filter(at->at.isCarbon())
    	.allMatch(c->c.isSp3Carbon()));
    }
    
    @Test
    public void testIsPiperazine() throws Exception {
    	FeaturizeNitrosamine.GLOBAL_SETTINGS.DO_EXTENDED_FEATURES_TOO=true;
    	
    	Chemical c1= Chemical.parse("N1CCNCC1");
    	    
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	assertEquals(2,resp.size());
    	FeatureResponse resp1= resp.get(0);
    	
    	assertEquals("A. Multiple Secondary Amine" ,resp1.getType());
    	assertEquals("2,2" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.ALPHA_HYDROGENS.getFeatureName()).orElse(null));
    	assertEquals("YES" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.PIPERAZINE.getFeatureName()).orElse(null));
    	
    	
    }
    
    @Test
    public void testCarboxylicAcidOnSaltDoesNotCount() throws Exception {
    	Chemical c1= Chemical.parse("O[C@H]([C@@H](O)C(O)=O)C(O)=O.COC1=CC=C(C[C@@H](C)[NH:20]C[C@H](O)C2=CC=C(O)C(NC=O)=C2)C=C1");
    	
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	assertEquals(1,resp.size());
    	FeatureResponse resp1= resp.get(0);
    	
    	assertEquals("A. Secondary Amine" ,resp1.getType());
    	assertEquals("NO" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.COOH.getFeatureName()).orElse(null));
//    	assertEquals("YES" ,resp1.getFeature("isPiperazine").orElse(null));
    	
    	
    }
    
    // DUVELISIB
    // IDELALISIB
    // OLUTASIDENIB
    // has potential benzyl, depending on aromatic model
    //
    
    @Test
    public void testBenzylLikeFeatureShouldNotFindPsuedoAromaticity() throws Exception {
    	Chemical c1= Chemical.parse("C[C@H]([NH:3]C1=C2N=CNC2=NC=N1)C3=CC4=C(C(Cl)=CC=C4)C(=O)N3C5=CC=CC=C5");
    	
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	assertEquals(1,resp.size());
    	FeatureResponse resp1= resp.get(0);
    	assertEquals("A. Secondary Amine" ,resp1.getType());
    	assertEquals("NO" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.ARYL_ALPHA.getFeatureName()).orElse(null));
    	
    }
    @Test
    public void testPseudoAromaticGroupStillCountsAsDoubleBond() throws Exception {

    	FeaturizeNitrosamine.GLOBAL_SETTINGS.DO_EXTENDED_FEATURES_TOO=true;
    	
    	Chemical c1= Chemical.parse("C[C@H]([NH:3]C1=C2N=CNC2=NC=N1)C3=CC4=C(C(Cl)=CC=C4)C(=O)N3C5=CC=CC=C5");
    	
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	assertEquals(1,resp.size());
    	FeatureResponse resp1= resp.get(0);
    	assertEquals("YES" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.ALLYL_GROUP.getFeatureName()).orElse(null));
    	
    }
    
    @Test
    public void testPseudoAromaticGroupStillCountsAsDoubleBondEWG() throws Exception {
    	Chemical c1= Chemical.parse("CN1C=C[NH:5]C1=S");
    	
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	assertEquals(1,resp.size());
    	FeatureResponse resp1= resp.get(0);
    	assertEquals("YES" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.EWG_ONE_SIDE.getFeatureName()).orElse(null));
    	
    }
    
    @Test
    public void testExocyclicAlleneNotAromatic() throws Exception {
    	Chemical c1= Chemical.parse("O=C1C=CNc2ccccc12");
    	
    	c1.aromatize();
    	Set<AtomTest> carbonylAromatic = AtomTest.stream(c1)
				.filter(at->at.isCarbonyl())
				.filter(cat->cat.isAromatic())
				.flatMap(at->at.getNeighbors(nat->nat.isOxygen()))
				.collect(Collectors.toSet());	
//    	long countBefore = BondTest.stream(c1).filter(b->b.isAromatic()).count();
//    	assertEquals(11,countBefore);
		carbonylAromatic.forEach(at->{
			at.get().setAtomicNumber(6);
		});
		c1.aromatize();
    	long countAfter = BondTest.stream(c1).filter(b->b.isAromatic()).count();
    	assertEquals(6,countAfter);
    	
    }
    
    @Test
    public void testFumerateSaltFind() throws Exception {
    	Chemical c1= Chemical.parse("OC(=O)\\C=C\\C(O)=O.COC1=NC2=CC=C(Br)C=C2C=C1[C@@H](C3=CC=CC=C3)[C@@](O)(CC[N:33](C)C)C4=C5C=CC=CC5=CC=C4");
    	
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	assertEquals(1,resp.size());
    	FeatureResponse resp1= resp.get(0);
    	assertEquals("Fumeric Acid", resp1.getSaltInfo().getName());
    	assertEquals("Salt", resp1.getSaltInfo().getIncludeType());
    	
    }
    
    
    @Test
    public void testLargeRingExcluded() throws Exception {
    	Chemical c1= Chemical.parse("N1CCCCCCCC1");
    	
    	long ringNs = AtomTest.stream(c1)
    	.filter(cat->cat.isNitrogen())
    	.filter(cat->cat.isInRing())
    	.count();
    	
    	assertEquals(1,ringNs);
    	
    }
    
    @Test
    public void testLayoutWorksOnCertainSalts() throws Exception {
    	Chemical c1= Chemical.parse("[Na+].O=C1N=C[NH]C2=C1C=N[N-]2");
	    
    	FeatureJob fj = new FeatureJob(c1);
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);

    	FeatureResponse resp1= resp.get(0);
    	assertEquals("Sodium", resp1.getSaltInfo().getName());
    	assertEquals("Salt", resp1.getSaltInfo().getIncludeType());

    }
    
    @Test
    public void testCOOHIsNoOnALLOPURINOL() throws Exception {
    	Chemical c1= Chemical.parse("O=c1nc[nH]c2c1cn[nH]2");
	    FeatureJob fj = new FeatureJob(c1);
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);    	
    	FeatureResponse resp1= resp.get(0);
    	assertEquals("NO" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.COOH.getFeatureName()).orElse(null));
    	

    }
    
    //


    @Test
    public void testZincSaltParsedOut() throws Exception {
    	Chemical c1= Chemical.parse("[Zn++].[H][C@](NC(=O)[C@@H](CCC([O-])=O)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H]1CSC(=N1)[C@@H](N)[C@@H](C)CC)([C@@H](C)CC)C(=O)N[C@H]2CCCCNC(=O)[C@H](CC(N)=O)NC(=O)[C@@H](CC([O-])=O)NC(=O)[C@H](CC3=C[NH:70]C=N3)NC(=O)[C@@H](CC4=CC=CC=C4)NC(=O)[C@@]([H])(NC(=O)[C@@]([H])(CCCN)NC2=O)[C@@H](C)CC");
    

    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	FeatureResponse resp1= resp.get(0);
    	
    	assertEquals(1,resp.size());
    	assertEquals("Zinc", resp1.getSaltInfo().getName());
    	assertEquals("Salt", resp1.getSaltInfo().getIncludeType());
    }
    
    @Test
    public void testNonAromaticDoubleBondHeteroExcluded() throws Exception {
    	Chemical c1= Chemical.parse("Nc1nc(=O)c2c(NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)=N2)[nH]1");
    	    
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	assertEquals(2,resp.size());
    }
    @Test
    public void testAbacavirSulfateSaltParsedOut() throws Exception {
    	Chemical c1= Chemical.parse("OS(O)(=O)=O.NC1=NC2=C(N=CN2[C@@H]3C[C@H](CO)C=C3)C(NC4CC4)=N1");
    	    
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	FeatureResponse resp1= resp.get(0);
    	
    	assertEquals(1,resp.size());
    	assertEquals("Sulfate", resp1.getSaltInfo().getName());
    	assertEquals("Salt", resp1.getSaltInfo().getIncludeType());
    }

    @Test
    public void testAbacavirSulfateSecondaryAmine() throws Exception {
    	Chemical c1= Chemical.parse("OS(O)(=O)=O.NC1=NC2=C(N=CN2[C@@H]3C[C@H](CO)C=C3)C(NC4CC4)=N1.NC5=NC6=C(N=CN6[C@@H]7C[C@H](CO)C=C7)C([NH:43]C8CC8)=N5");
    	    
    	FeatureJob fj = new FeatureJob(c1);
    	
    	List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);
    	
    	FeatureResponse resp1= resp.get(0);
    	
    	assertEquals(1,resp.size());
    	assertEquals("A. Secondary Amine" ,resp1.getType());
    	assertEquals("0,1" ,resp1.getFeature(FeaturizeNitrosamine.FeaturePairRegistry.ALPHA_HYDROGENS.getFeatureName()).orElse(null));
    }
	//NC1=NC2=C(N=CN2[C@@H]3C[C@H](CO)C=C3)C(NC4CC4)=N1

	@Test
	public void testTartaricAcid() throws Exception {
		//we expect this molecule not to have any nitrosamine potential whatsoever
		Chemical c1= Chemical.parse("[C@@H]([C@H](C(=O)O)O)(C(=O)O)O");

		FeatureJob fj = new FeatureJob(c1);

		List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);

		assertNotNull(resp);
		System.out.println("got output");
		assertTrue( resp.isEmpty());
	}

	@Test
	public void testAmide() throws Exception {
		//we expect this molecule not to have any nitrosamine potential whatsoever
		Chemical c1= Chemical.parse("C(C)CCN(C(N)=N)N=O");

		FeatureJob fj = new FeatureJob(c1);

		List<FeatureResponse> resp = FeaturizeNitrosamine.fingerprintNitrosamine(fj);

		assertNotNull(resp);
		System.out.println("got output");
		assertTrue( resp.isEmpty());
	}

	@Test
	public void testNitrosamine() throws Exception {
		Chemical testChemical = Chemical.parse("CCN(CC)N=O");
		FeatureJob job = new FeatureJob(testChemical);
		List<FeatureResponse> responses = FeaturizeNitrosamine.fingerprintNitrosamine(job);
		assertNotNull(responses);
		for(FeatureResponse resp: responses) {
			System.out.printf("type: %s\n", resp.getType());
		}
	}

	@Test
	public void testNitrosamideSetFromMolfile() throws Exception {
		List<TestMol> testData = Arrays.asList(new TestMol("nitrosamiide-like", false,
				"nitrosamide-like with S expects no output"),
				new TestMol("nitrosamiide-like2", false,
						"nitrosamide-like with N expects no output"),
				new TestMol("nitrosamine_amide", false, "amides are currently not supported by model"),
				new TestMol("nitrosamine1", true, "nitrosamine produces output"));
		testData.forEach(m->{
			try {
				InputStream stream = getClass().getResourceAsStream("/molfiles/" + m.molFileName + ".mol");
				String molfileText = IOUtils.toString(stream, "UTF-8");
				Chemical testChemical = Chemical.parse(molfileText);
				Optional<FeatureResponse> response = FeaturizeNitrosamine.forMostPotentNitrosamine(testChemical);
				assertEquals(m.expectData, response.isPresent());
			} catch (Exception ex){
				Logger.getLogger(this.getClass().getName()).warning(String.format("failure for %s", m.description));
				fail();
			}
		});
		/*InputStream stream =getClass().getResourceAsStream("/molfiles/nitrosamiide-like.mol");
		String molfileText = IOUtils.toString(stream, "UTF-8");
		Chemical testChemical = Chemical.parse(molfileText);
		Optional<FeatureResponse> response = FeaturizeNitrosamine.forMostPotentNitrosamine(testChemical);
		assertFalse(response.isPresent());*/
	}

}
