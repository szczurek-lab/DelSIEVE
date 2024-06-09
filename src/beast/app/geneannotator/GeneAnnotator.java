package beast.app.geneannotator;

import beast.app.tools.LogCombiner;
import beast.app.util.Arguments;
import beast.app.util.Utils;
import beast.app.utils.ChromosomeLabel;
import beast.core.util.Log;
import beast.evolution.alignment.VariantSiteInfo;
import beast.evolution.substitutionmodel.ScsFiniteMuDelModel;
import beast.evolution.substitutionmodel.ScsFiniteMuExtendedModel;
import beast.evolution.substitutionmodel.ScsSubstitutionModelBase;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.NexusParser;
import beast.util.NotSingleException;
import beast.util.TreeParser;
import jam.console.ConsoleApplication;
import org.jetbrains.annotations.NotNull;

import javax.swing.*;
import java.io.*;
import java.nio.file.Files;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static beast.evolution.variantsinfo.GenericVariantsInfoVCF.META_DATA_GENOTYPES;
import static beast.util.FileNameProcessor.getBaseName;
import static beast.util.FileNameProcessor.getRootPath;
import static beast.util.TreeUtils.processMetaData;

public class GeneAnnotator {

    private final static String GENOTYPES_FORMAT_ERROR = "The genotypes parsed from the input tree can only be numbers.";

    private final static String METADATA_INDEX = "index";
    private final static String METADATA_CHR = "chr";
    private final static String METADATA_POS = "pos";
    private final static String METADATA_REF_NUC = "ref_nuc";
    private final static String METADATA_ALT_NUC = "alt_nuc";
    private final static String METADATA_EVENT_TYPE = "event_type";
    private final static String METADATA_GENE = "gene";
    private final static String METADATA_ISA_SUFFIX = "_isa";
    private final static String METADATA_FSA_SUFFIX = "_fsa";

    private final static String OUTPUT_SUFFIX = ".gene_tree";

    private final static String TREE_PATTERN = ".+\\.intermediate_tree$";

    private final static String MUTATED_SITE_PATTERN = ".+\\.loci_info$";

    private final static String GT_TERNARY_PATTERN = ".+\\.ternary$";

    private final static String GENE_NAME_HEADER = "Gene.refGene";

    static PrintStream progressStream = Log.err;

    private NodeGenesInfo[] nodeGenesInfos;


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public GeneAnnotator(
            int substModelLabel,
            String path2results,
            String mutationMapFileName,
            boolean mapFromAnnovar,
            String mapSeparator,
            int annovarFilterCol,
            String[] annovarFilterKeywords,
            String filteringGenesFileName,
            String filteringGenesSeparator,
            int filteringGenesColIndex,
            String outputFileName
    ) throws IOException, NotSingleException {
        ScsSubstitutionModelBase substModel;
        final File path = new File(path2results);
        Tree tree;
        List<VariantSiteInfo> snvSites;
        List<Boolean> hasCopyNumberChanges = null;
        List<String> filteringGenes = null;
        Map<ChromosomeLabel, List<MutationMap>> mutationMap;

        // Sanity check.
        if (!path.isDirectory()) {
            Log.err.println(path2results + " does not exist!");
            return;
        }

        // 1. Get substitution model
        substModel = getSubstModel(substModelLabel);

        // 2. Get input tree
        tree = getTree(
                getFile(
                        path,
                        TREE_PATTERN,
                        1
                )
        );

        // 3. Get variant sites info
        snvSites = getSNVSitesInfo(
                getFile(
                        path,
                        MUTATED_SITE_PATTERN,
                        1
                )
        );

        // 4. Get filtering genes
        if (filteringGenesFileName != null)
            filteringGenes = getFilteringGenes(filteringGenesFileName, filteringGenesSeparator, filteringGenesColIndex);

        // 5. Get mutation map
        if (mapFromAnnovar) {
            List<VariantSiteInfo> sitesToKeep = null;
            if (annovarFilterCol >= 0)
                sitesToKeep = getSitesToKeep(
                        getFile(
                                path,
                                GT_TERNARY_PATTERN,
                                1
                        ),
                        snvSites
                );

            mutationMap = getMutationMapFromAnnovar(
                    mutationMapFileName,
                    mapSeparator,
                    annovarFilterCol,
                    annovarFilterKeywords,
                    sitesToKeep,
                    filteringGenes
            );
        } else
            mutationMap = getMutationMap(mutationMapFileName, filteringGenes);

        // 6. Initialize some variables
        this.nodeGenesInfos = new NodeGenesInfo[tree.getNodeCount()];

        // 7. Collect genes for all nodes
        collectGenes4Tree(
                tree.getRoot(),
                META_DATA_GENOTYPES,
                null,
                substModel,
                snvSites,
                mutationMap
        );

        // 8. Process genes for all nodes
        processGenes();

        // 9. Annotate the tree
        for (Node node : tree.getNodesAsArray()) {
            node.removeMetaData(META_DATA_GENOTYPES);

            final NodeGenesInfo info = nodeGenesInfos[node.getNr()];

            if (info == null) continue;

            node.setMetaData(METADATA_INDEX, node.getNr());

            if (info.getFullGenes() != null && info.getFullGenes().size() > 0) {
                node.setMetaData(METADATA_CHR, info.getFullChr().toArray(new String[0]));
                node.setMetaData(METADATA_POS, info.getFullPos().toArray(new Long[0]));
                node.setMetaData(METADATA_REF_NUC, info.getFullRefNuc().toArray(new String[0]));
                node.setMetaData(METADATA_ALT_NUC, info.getFullAltNuc().toArray(new String[0]));
                node.setMetaData(METADATA_EVENT_TYPE, info.getFullEventType().toArray(new String[0]));
                node.setMetaData(METADATA_GENE, info.getFullGeneName().toArray(new String[0]));
            }

            if (info.getISAGenes() != null && info.getISAGenes().size() > 0) {
                node.setMetaData(METADATA_CHR + METADATA_ISA_SUFFIX, info.getISAChr().toArray(new String[0]));
                node.setMetaData(METADATA_POS + METADATA_ISA_SUFFIX, info.getISAPos().toArray(new Long[0]));
                node.setMetaData(METADATA_REF_NUC + METADATA_ISA_SUFFIX, info.getISARefNuc().toArray(new String[0]));
                node.setMetaData(METADATA_ALT_NUC + METADATA_ISA_SUFFIX, info.getISAAltNuc().toArray(new String[0]));
                node.setMetaData(METADATA_EVENT_TYPE + METADATA_ISA_SUFFIX, info.getISAEventType().toArray(new String[0]));
                node.setMetaData(METADATA_GENE + METADATA_ISA_SUFFIX, info.getISAGeneName().toArray(new String[0]));
            }

            if (info.getFSAGenes() != null && info.getFSAGenes().size() > 0) {
                node.setMetaData(METADATA_CHR + METADATA_FSA_SUFFIX, info.getFSAChr().toArray(new String[0]));
                node.setMetaData(METADATA_POS + METADATA_FSA_SUFFIX, info.getFSAPos().toArray(new Long[0]));
                node.setMetaData(METADATA_REF_NUC + METADATA_FSA_SUFFIX, info.getFSARefNuc().toArray(new String[0]));
                node.setMetaData(METADATA_ALT_NUC + METADATA_FSA_SUFFIX, info.getFSAAltNuc().toArray(new String[0]));
                node.setMetaData(METADATA_EVENT_TYPE + METADATA_FSA_SUFFIX, info.getFSAEventType().toArray(new String[0]));
                node.setMetaData(METADATA_GENE + METADATA_FSA_SUFFIX, info.getFSAGeneName().toArray(new String[0]));
            }
        }

        // 10. Save the annotated tree to disk
        processMetaData(tree.getRoot());
        final PrintStream treeOut = new PrintStream(outputFileName);
        tree.init(treeOut);
        treeOut.println();
        treeOut.print("tree TREE1 = ");
        treeOut.print(tree.getRoot().toSortedNewick(new int[1], true));
        treeOut.println(";");
        tree.close(treeOut);
        treeOut.println();
        treeOut.close();

        progressStream.println("Successful!");
        progressStream.println("Annotated tree has been save to " + outputFileName);
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    private ScsSubstitutionModelBase getSubstModel(final int substModelLabel) {
        switch (substModelLabel) {
            case 0:
                return new ScsFiniteMuExtendedModel();
            case 1:
                return new ScsFiniteMuDelModel();
            default:
                throw new RuntimeException("Error! Invalid label for substitution model: " + substModelLabel);
        }
    } // getSubstModel

    private String getFile(@NotNull final File path, final String pattern, final int num_expected) {
        final Pattern pat = Pattern.compile(pattern);

        final String[] children = path.list();
        assert children != null: "Error! The folder is empty: " + path;

        final List<String> founds = Arrays.stream(children).filter(x -> pat.matcher(x).matches()).collect(Collectors.toList());
        assert founds.size() == num_expected: "Error! There should be " + num_expected + " existing file in " + path + " matching this pattern " + pattern;

        return new File(path, founds.get(0)).getPath();
    } // getFile

    private Tree getTree(final String treeFileName) throws IOException, NotSingleException {
        List<Tree> parsedTrees;
        boolean isNexus = true;
        BufferedReader fin;
        String str;

        // Nexus tree or Newick tree?
        fin = new BufferedReader(new FileReader(treeFileName));
        if (!fin.ready()) {
            throw new IOException(treeFileName + " appears empty.");
        }
        str = fin.readLine();
        if (!str.toUpperCase().trim().startsWith("#NEXUS")) {
            // Newick tree found
            isNexus = false;
        }

        // Parse tree
        if (isNexus) {
            NexusParser nexusParser = new NexusParser();
            nexusParser.parseFile(new File(treeFileName));
            parsedTrees = nexusParser.trees;
        } else {
            parsedTrees = new ArrayList<>();
            while (fin.ready()) {
                str = fin.readLine().trim();

                Tree thisTree;
                try {
                    thisTree = new TreeParser(null, str, 0, false);
                } catch (ArrayIndexOutOfBoundsException e) {
                    thisTree = new TreeParser(null, str, 1, false);
                }

                parsedTrees.add(thisTree);
            }
        }
        fin.close();

        // Verify; only one input tree is allowed
        if (parsedTrees.size() != 1) {
            throw new NotSingleException("Only one input tree is expected, but " + parsedTrees.size() + " provided.");
        }

        return parsedTrees.get(0);
    } // getTree

    private List<VariantSiteInfo> getSNVSitesInfo(final String snvFileName) throws IOException {
        List<VariantSiteInfo> results = new ArrayList<>();
        BufferedReader fin = new BufferedReader(new FileReader(snvFileName));

        if (!fin.ready()) {
            throw new IOException(snvFileName + " appears empty.");
        }

        while (fin.ready()) {
            String[] comp = fin.readLine().trim().split("\\s+");

            try {
                results.add(
                        new VariantSiteInfo(
                                comp[0],
                                Long.parseLong(comp[1]),
                                comp[2],
                                comp[3].split(",")
                        )
                );
            } catch (ArithmeticException e) {
                throw new IllegalArgumentException("Error! Make sure the position is a number: " +
                        String.join(";", comp));
            }
        }

        fin.close();
        return results;
    } // getSNVSitesInfo

    /**
     * Get the variant sites to keep. For now, any sites containing copy number changes are kept.
     *
     * @param gtTernaryFileName apparently
     * @param snvSites          snv sites corresponding to the ternary genotypes
     * @return snv sites to keep.
     * @throws IOException
     */
    private List<VariantSiteInfo> getSitesToKeep(
            final String gtTernaryFileName,
            final List<VariantSiteInfo> snvSites
    ) throws IOException {
        assert new File(gtTernaryFileName).exists();
        List<int[]> gts = new ArrayList<>();
        List<VariantSiteInfo> results = new ArrayList<>();

        BufferedReader fin = new BufferedReader(new FileReader(gtTernaryFileName));
        if (!fin.ready()) {
            throw new IOException(gtTernaryFileName + " appears empty.");
        }

        while (fin.ready()) {
            String[] comp = fin.readLine().trim().split("\\s+");

            try {
                gts.add(
                        Arrays.stream(comp).mapToInt(Integer::parseInt).toArray()
                );
            } catch (ArithmeticException e) {
                throw new IllegalArgumentException("Error! Make sure the position is a number: " +
                        String.join(";", comp));
            }
        }
        fin.close();

        assert gts.size() == snvSites.size(): "Error! The number of lines in the SNV sites and ternary genotypes should be the same.";

        for (int i = 0; i < gts.size(); i++) {
            final int[] hasCNA = Arrays.stream(gts.get(i)).filter(x -> x < 0 || x > 3).toArray();

            if (hasCNA.length > 0)
                results.add(snvSites.get(i));
        }

        return results;
    } // getSitesToKeep

    private List<String> getFilteringGenes(
            String filteringGenesFileName,
            String filteringGenesSeparator,
            int filteringGenesColIndex
    ) throws IOException {
        if (filteringGenesFileName == null)
            return null;

        List<String> results = new ArrayList<>();
        BufferedReader fin = new BufferedReader(new FileReader(filteringGenesFileName));

        if (!fin.ready())
            throw new IOException(filteringGenesFileName + " appears empty.");

        int lineNumber = 0;
        while (fin.ready()) {
            final String[] line = fin.readLine().trim().replaceAll("\"", "").split(filteringGenesSeparator);

            if (lineNumber > 0) {
                final String gene = line[filteringGenesColIndex];

                if (!results.contains(gene))
                    results.add(gene);
            }

            lineNumber++;
        }

        fin.close();

        return results;
    } // getFilteringGenes

    private Map<ChromosomeLabel, List<MutationMap>> getMutationMapFromAnnovar(
            final String mutationMapFileName,
            final String mapSeparator,
            final int annovarFilterCol,
            final String[] annovarFilterKeywords,
            final List<VariantSiteInfo> sitesToKeep,
            final List<String> filteringGenes
    ) throws IOException {
        Map<ChromosomeLabel, List<MutationMap>> results = new HashMap<>();
        assert new File(mutationMapFileName).isFile();
        BufferedReader fin = new BufferedReader(new FileReader(mutationMapFileName));

        if (!fin.ready())
            throw new IOException(mutationMapFileName + " appears empty.");

        int lineNumber = 0;
        int geneNameIndex = -1;
        while (fin.ready()) {
            final String[] line = fin.readLine().trim().replaceAll("\"", "").split(mapSeparator);

            if (lineNumber == 0) {
                // process the header

                for (int i = 0; i < line.length; i++) {
                    if (line[i].equalsIgnoreCase(GENE_NAME_HEADER)) {
                        geneNameIndex = i;
                        break;
                    }
                }

                assert annovarFilterCol < line.length: "Error! The column index used to filter genes in the results of Annovar should not exceed " + line.length;

                if (geneNameIndex == -1)
                    throw new IOException("Error! The mutation map does not contain a column named " + GENE_NAME_HEADER);
            } else {
                // process each line below the header

                MutationMap mm = new MutationMap(new String[]{line[0].trim(), line[1].trim(), line[2].trim(), line[geneNameIndex].trim()});
                final ChromosomeLabel key = mm.getChromosome();

                if (
                        (filteringGenes == null || filteringGenes.contains(mm.getGeneName())) &&
                                (annovarFilterCol < 0 || keepGene(line[annovarFilterCol], annovarFilterKeywords) ||
                                        sitesToKeep == null || sitesToKeep.stream().anyMatch(mm::containsVariantSite))
                ) {
                    if (results.containsKey(key))
                        results.get(key).add(mm);
                    else
                        results.put(key, Stream.of(mm).collect(Collectors.toList()));
                }
            }

            lineNumber++;
        }

        final MutationMap.ListComparator comparator = new MutationMap.ListComparator();
        for (ChromosomeLabel chromosomeLabel : results.keySet()) {
            results.get(chromosomeLabel).sort(comparator);
        }

        fin.close();
        return results;
    } // getMutationMapFromAnnovar

    private boolean keepGene(final String test, final String[] keywords) {
        assert keywords.length > 1: "Error! Option -anvfilterkws has a wrong format.";

        if (Objects.equals(keywords[0], "r")) {
            final Pattern pat = Pattern.compile(
                    String.join("|", Arrays.copyOfRange(keywords, 1, keywords.length))
            );

            return pat.matcher(test).matches();
        } else if (Objects.equals(keywords[0], "f")) {
            for (int i = 1; i < keywords.length; i++)
                if (Objects.equals(test, keywords[i]))
                    return true;

            return false;
        }

        throw new RuntimeException("Error! Option -anvfilterkws should start with either 'r' or 'f'.");
    } // keepGene

    private Map<ChromosomeLabel, List<MutationMap>> getMutationMap(
            String mutationMapFileName,
            List<String> filteringGenes
    ) throws IOException {
        Map<ChromosomeLabel, List<MutationMap>> results = new HashMap<>();
        BufferedReader fin = new BufferedReader(new FileReader(mutationMapFileName));

        if (!fin.ready())
            throw new IOException(mutationMapFileName + " appears empty.");

        while (fin.ready()) {
            String line = fin.readLine().trim();
            if (!line.startsWith("#")) {
                MutationMap mm = new MutationMap(line.split("\\s+"));
                final ChromosomeLabel key = mm.getChromosome();

                if (filteringGenes == null || filteringGenes.contains(mm.getGeneName())) {
                    if (results.containsKey(key))
                        results.get(key).add(mm);
                    else
                        results.put(key, Stream.of(mm).collect(Collectors.toList()));
                }
            }
        }

        final MutationMap.ListComparator comparator = new MutationMap.ListComparator();
        for (ChromosomeLabel chromosomeLabel : results.keySet()) {
            results.get(chromosomeLabel).sort(comparator);
        }

        fin.close();
        return results;
    } // getMutationMap

    private void collectGenes4Tree(
            Node node,
            String genotypeKey,
            Object parentGenotypes,
            final ScsSubstitutionModelBase substModel,
            final List<VariantSiteInfo> snvSites,
            final Map<ChromosomeLabel, List<MutationMap>> mutationMap
    ) throws IllegalArgumentException {
        if (!node.isRoot() && !(parentGenotypes instanceof Integer[]) && !(parentGenotypes instanceof Double[]))
            throw new IllegalArgumentException(GENOTYPES_FORMAT_ERROR);

        Object genotypes = node.getMetaData(genotypeKey);
        if (node.isLeaf() && !(genotypes instanceof Integer[]) && !(genotypes instanceof Double[]))
            throw new IllegalArgumentException(GENOTYPES_FORMAT_ERROR);

        // Walk the tree in preorder
        for (Node child : node.getChildren())
            collectGenes4Tree(child, genotypeKey, genotypes, substModel, snvSites, mutationMap);

        if (!node.isRoot()) {
            for (int snv = 0; snv < snvSites.size(); snv++) {
                String geneName = getGeneName(snvSites.get(snv), mutationMap);

                // This snv site is located in the range of an important gene.
                if (geneName != null) {
                    final int childGenotype = (((Double[]) genotypes)[snv]).intValue();
                    final int parentGenotype = (((Double[]) parentGenotypes)[snv]).intValue();

                    if (parentGenotype != childGenotype) {
                        if (nodeGenesInfos[node.getNr()] == null)
                            nodeGenesInfos[node.getNr()] = new NodeGenesInfo();

                        ScsSubstitutionModelBase.EvolutionaryEventType[] events = substModel.getEvolutionaryEvents(
                                parentGenotype,
                                childGenotype
                        );

                        if (events != null && events.length > 0) {
                            // replace ';' with '/' in geneName to avoid subsequent parsing problems
                            geneName = geneName.replaceAll(";", "/");

                            nodeGenesInfos[node.getNr()].addGeneInfo(
                                    new GeneInfo(
                                            snvSites.get(snv),
                                            events,
                                            geneName
                                    )
                            );
                        }
                    }
                }
            }
        }
    } // annotateGenes2Tree

    private String getGeneName(final VariantSiteInfo info, final Map<ChromosomeLabel, List<MutationMap>> mutationMap) {
        final ChromosomeLabel cl = info.getChromosomeLabel();
        if (mutationMap.containsKey(cl)) {
            for (MutationMap map : mutationMap.get(cl)) {
                final long pos = info.getPosition();
                if (pos >= map.getStartPos() && pos <= map.getEndPos())
                    return map.getGeneName();
            }
        }

        return null;
    } // getGeneName

    private void processGenes() {
        List<GeneInfo> genes = new ArrayList<>();
        List<GeneInfo> genesISA = new ArrayList<>();
        List<GeneInfo> genesFSA = new ArrayList<>();

        for (NodeGenesInfo i : this.nodeGenesInfos) {
            if (i != null && i.getFullGenes() != null)
                genes.addAll(i.getFullGenes());
        }

        genes.sort(new GeneInfo.GeneInfoComparator());

        for (GeneInfo i : new HashSet<>(genes)) {
            final int j = Collections.frequency(genes, i);
            if (j > 1 || (j == 1 && i.violateISA()))
                genesFSA.add(i);
            else
                genesISA.add(i);
        }

        for (NodeGenesInfo i : this.nodeGenesInfos) {
            if (i != null && i.getFullGenes() != null) {
                for (GeneInfo j : i.getFullGenes()) {
                    if (genesISA.contains(j)) {
                        i.addISAGeneInfo(j);
                    } else if (genesFSA.contains(j)) {
                        i.addFSAGeneInfo(j);
                    }
                }

                i.convertGeneInfo();
            }
        }
    } // processGenes


    //**********************************************
    //*               Static methods               *
    //**********************************************

    private static String getOutputFileName(final String fileName) {
        return getRootPath(fileName) + getBaseName(fileName) + OUTPUT_SUFFIX;
    } // getOutputFileName

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("geneannotator", "");
    } // printUsage


    //***********************************************
    //*                 Main method                 *
    //***********************************************

    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        // Variables
        int substModelLabel;
        String path2results;
        String mutationMapFileName;
        String outputFileName;
        boolean mapFromAnnovar = false;
        int annovarFilterCol = -1;
        String[] annovarFilterKeywords = null;
        int mapSeparator = -1;
        String filteringGenesFileName = null;
        int filteringGenesSeparator = 0;
        int filteringGenesColIndex = 0;

        if (args.length == 0) {

            Utils.loadUIManager();
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");
            java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
            javax.swing.Icon icon = null;
            if (url != null)
                icon = new javax.swing.ImageIcon(url);

            // Construct a new console
            new ConsoleApplication(null, null, icon, true);
            Log.info = System.out;
            Log.err = System.err;
            progressStream = System.out;

            // TODO: print some information here
            System.out.println("GeneAnnotator");

            GeneAnnotatorDialog dialog = new GeneAnnotatorDialog(new JFrame());
            if (!dialog.showDialog("GeneAnnotator"))
                return;

            // Get parameters
            path2results = dialog.getPath2results();
            if (path2results == null) {
                Log.err.println("No path to the results of VariantCaller specified!");
                return;
            }

            substModelLabel = dialog.getSubstModelLabel();

            mutationMapFileName = dialog.getMutationMapFileName();
            if (mutationMapFileName == null) {
                Log.err.println("No mutation map file specified!");
                return;
            }

            mapFromAnnovar = dialog.isMapFromAnnovar();
            if (mapFromAnnovar)
                mapSeparator = dialog.getMapSeparator();

            annovarFilterCol = dialog.getAnnovarFilterColIndex();
            annovarFilterKeywords = dialog.getAnnovarFilteringKeywords();

            filteringGenesFileName = dialog.getFilteringFileName();
            filteringGenesSeparator = dialog.getFilteringFileSeparator();
            filteringGenesColIndex = dialog.getColIndexFilteringGenes();

            outputFileName = dialog.getOutputFileName();
            if (outputFileName == null)
                outputFileName = getOutputFileName(path2results);

            try {
                new GeneAnnotator(
                        substModelLabel,
                        path2results,
                        mutationMapFileName,
                        mapFromAnnovar,
                        mapSeparator != 0 && mapSeparator != 1 ? null : (mapSeparator == 0 ? "\t" : ","),
                        annovarFilterCol,
                        annovarFilterKeywords,
                        filteringGenesFileName,
                        filteringGenesSeparator != 0 && filteringGenesSeparator != 1 ? null : (filteringGenesSeparator == 0 ? "\t" : ","),
                        filteringGenesColIndex, outputFileName);
            } catch (Exception e) {
                e.printStackTrace();
                progressStream.println("Exception: " + e.getMessage());
            }

            progressStream.println("Finished - Quit program to exit.");
            while (true) {
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

        } else {

            // TODO: print some information here
            Arguments arguments = new Arguments(
                    new Arguments.Option[]{
                            new Arguments.Option("help", "option to print this message -> OPTIONAL"),
                            new Arguments.StringOption("details", "variant_caller_details", "specifies a directory containing the output of VariantCaller called with option '-details'. -> MANDATORY"),
                            new Arguments.IntegerOption("subst", 0, 1, "specifies the substitution model used to infer phylogeny and call variants (0: ScsFiniteMuExtendedModel, 1: ScsFiniteMuDelModel). -> MANDATORY"),
                            new Arguments.StringOption("map", "mutation_map_file", "specifies the mutation map. -> MANDATORY"),
                            new Arguments.IntegerOption("anv", 0, 1, "specifies when the value of --map option is an output of Annovar (0: white-space-separated, 1: comma-separated). -> OPTIONAL"),
                            new Arguments.IntegerOption("anvfiltercol", -1, 50000, "specifies which column of the output of Annovar should be used to filter genes; only working when --anv is specified; setting to -1 for non-filtering (0-based; default: -1). -> OPTIONAL"),
                            new Arguments.StringOption("anvfilterkws", "annovar_filtering_keywords", "specifies the comma-separated keywords to keep genes, with the first keyword being either 'f' for 'fixed' or 'r' for 'regular expression'. -> OPTIONAL"),
                            new Arguments.StringOption("filter", "filtering_genes_file", "specifies filtering genes file (headers required). -> OPTIONAL"),
                            new Arguments.IntegerOption("sep", 0, 1, "specifies the separator of filtering genes file (0: tab-separated; 1: comma-separated; default: 0). -> OPTIONAL"),
                            new Arguments.IntegerOption("col", 0, 1000, "specifies the column index of gene names in the filtering genes file (0-based; default: 0). -> OPTIONAL"),
                            new Arguments.StringOption("out", "out_tree_file", "specifies the output tree file with genes annotated. -> OPTIONAL")
                    }
            );

            // Parse command line arguments
            try {
                arguments.parseArguments(args);
            } catch (Arguments.ArgumentException e) {
                progressStream.println(e);
                printUsage(arguments);
                System.exit(1);
            }

            // Print help message
            if (arguments.hasOption("help") || arguments.hasOption("h")) {
                printUsage(arguments);
                System.exit(0);
            }

            // Set path2results
            if (arguments.hasOption("details"))
                path2results = arguments.getStringOption("details");
            else {
                Log.err.println("No path to the results of VariantCaller specified!");
                return;
            }

            // Set substModelLabel
            if (arguments.hasOption("subst"))
                substModelLabel = arguments.getIntegerOption("subst");
            else {
                Log.err.println("No substitution model specified! Check the usage.");
                return;
            }

            // Set mutationMapFileName
            if (arguments.hasOption("map"))
                mutationMapFileName = arguments.getStringOption("map");
            else {
                Log.err.println("No mutation map file specified!");
                return;
            }

            // Set mapFromAnnovar and mapSeparator
            if (arguments.hasOption("anv")) {
                mapFromAnnovar = true;
                mapSeparator = arguments.getIntegerOption("anv");
            }

            // Set annovarFilterCol
            if (arguments.hasOption("anvfiltercol")) {
                annovarFilterCol = arguments.getIntegerOption("anvfiltercol");
            }

            // Set annovarFilterKeywords
            if (arguments.hasOption("anvfilterkws")) {
                annovarFilterKeywords = arguments.getStringOption("anvfilterkws").split(",");
            }

            // Sanity check of Annovar-related options
            if ((annovarFilterCol >= 0 || annovarFilterKeywords != null) && !mapFromAnnovar) {
                Log.err.println("'--anv' must be specified to have '--anvfiltercol' and '--anvfilterkws' come into effect.");
                return;
            } else if (annovarFilterCol >= 0 && annovarFilterKeywords == null || annovarFilterCol < 0 && annovarFilterKeywords != null) {
                Log.err.println("'--anvfiltercol' and '--anvfilterkws' must be specified at the same time.");
                return;
            }

            // Set filteringGenesFileName
            if (arguments.hasOption("filter"))
                filteringGenesFileName = arguments.getStringOption("filter");

            // Set filteringGenesSeparator
            if (arguments.hasOption("sep"))
                filteringGenesSeparator = arguments.getIntegerOption("sep");

            // Set filteringGenesColIndex
            if (arguments.hasOption("col"))
                filteringGenesColIndex = arguments.getIntegerOption("col");

            // Set outputFileName
            if (arguments.hasOption("out"))
                outputFileName = arguments.getStringOption("out");
            else
                outputFileName = getOutputFileName(path2results);

            try {
                new GeneAnnotator(
                        substModelLabel,
                        path2results,
                        mutationMapFileName,
                        mapFromAnnovar,
                        mapSeparator != 0 && mapSeparator != 1 ? null : (mapSeparator == 0 ? "\t" : ","),
                        annovarFilterCol,
                        annovarFilterKeywords,
                        filteringGenesFileName,
                        filteringGenesSeparator != 0 && filteringGenesSeparator != 1 ? null : (filteringGenesSeparator == 0 ? "\t" : ","),
                        filteringGenesColIndex, outputFileName);
            } catch (IOException e) {
                throw e;
            } catch (Exception e) {
                e.printStackTrace();
            }

            progressStream.println();
            progressStream.println("Successful!");
        }

    } // main

}
