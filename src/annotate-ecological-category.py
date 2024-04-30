#!/usr/bin/env python

'''
annotate-ecological-category.py by Vincent Huang and Rohan Maddamsetti.

This script annotates the ecology for each analyzed genome, based on the metadata in the
'host' and 'isolation_source' tags in its genome file from Genbank.

Some definitions: soil is everything from dirt, whether an agricultural field or not.
agriculture, on the other hand, refers strictly to crop plants-- soil is not included.

In addition, "rhizosphere" always maps to the Soil annotation.

The "Animal" category also includes fungi and protists-- basically all non-plant eukaryotes.

Usage: python annotate-ecological-category.py > ../results/computationally-annotated-gbk-annotation-table.csv
'''

#############################################################################################
## Keywords for Annotating Ecological Categories.
#############################################################################################

animal_hosts = [
    "hoplopleura acanthopus", "capra aegagrus", "iguana iguana", "phasianus colchicus",
    "eurypanopeus depressus", "dermanyssus gallinae (poultry red mite)", "polytelis", "jguana",
    "reticulitermes santonensis (termite)", "steinernema", "mus musculus c3h/hen", "monkey",
    "haemaphysalis megaspinosa", "heterorhabditis brevicaudis", "pineus similis", "ground squirrel",
    "rhagoletis cingulata", "ambystoma andersoni", "nysius blackburni", "nysius rubescens",
    "adelges cooleyi on picea engelmanii", "cacospongia mycofijiensis", "ixodes pavlovskyi",
    "scarabaeidae_pupa_coleoptera", "ailurus fulgens", "neovison vison", "spodoptera picta",
    "sturnus vulgaris", "amblyomma nuttalli", "wiebesia pumilae", "lamellibrachia satsuma 1372_tw2",
    "ixodes nipponensis", "sciurus vulgaris", "ailuropoda melanoleuca", "mosquito", "hamster",
    "philander opossum", "dermacentor parumapertus", "philanthus triangulum", "bats", "termite",
    "naemorhedus caudatus", "lithobius_myriapod", "thymallus thymallus", "cimex lectularius",
    "nysius terrestris", "elephas maximus", "seastar", "aedes aegypti", "bullfrog", "pagrus major",
    "drosophila simulans madagaskar", "dreyfusia piceae", "rattus rattus", "empoasca decipiens",
    "biomphalaria glabrata", "heosemys grandis", "starfish", "drosophila simulans stc", "stenella coeruleoalba",
    "lipoptena fortisetosa", "tursiops truncatus", "magicicada tredecim", "cacatua galerita",
    "thais luteostoma", "chondrilla nucula", "cassida rubiginosa", "shelfordina sp.", "lophura swinhoii",
    "epilampra maya", "xiphinema americanum", "rhipicephalus turanicus", "dendrocephalus brasiliensis",
    "cruorifilaria tuberocauda", "pyractomena angulata", "opisthoplatia orientalis", "tortoise",
    "mus musculus c57bl/6", "yezoterpnosia nigricosta", "anopheles gambiae lab colony",
    "batillaria multiformis", "elk", "crucian", "llama", "cryptotympana facialis", "bird (pigeon)",
    "camponotus pennsylvanicus", "therioaphis trifolii", "harpegnathos saltator", "onchocerca ochengi",
    "linaria flavirostris", "lolium rigidum and anguina funesta", "liza haematocheila", "adelges piceae",
    "african bull frog; pyxicephalus edulis", "trachymyrmex cornetzi", "meimuna iwasakii",
    "spheniscus demersus", "dsm1974", "graptopsaltria nigrofuscata", "castor sp.", "the pantholops hodgsonii",
    "perca fluviatilis", "vespid wasp", "muda kuroiwae", "tibetan antilope", "methana sp.",
    "acanthamoeba sp. atcc pra2", "camponotus floridanus", "rhopalosiphum maidis", "macaca fuscata",
    "omphisa fuscidentalis hampson", "polyrhachis (hedomyrma) turneri", "silurus asotus",
    "auritibicen japonicus", "stephanolepis cirrhifer", "cinara tujafilina", "tettigades auropilosa",
    "seriola quinqueradiata", "uroleucon sonchi", "angomonas deanei atcc 30255", "dalbulus maidis",
    "tanna japonensis", "sparus aurata", "haematobia irritans", "mice (c57bl/6)", "mammal (ferret)",
    "hartmannella sp. fs5", "puffer fish", "hoplobatrachus chinensis", "aegypius", "artipe eryx",
    "brugia pahangi", "sturnus nigricollis", "atrichopogon sp.", "lates calcarifer", "rupornis magnirostris",
    "hoplobatrachus rugulosus", "neogale vison", "chipmunk", "blackfly", "lizard (cordylus niger)",
    "dictyophara multireticulata", "carp", "pyrrhocoris apterus (red soldier bug)", "tachysurus fulvidraco",
    "paramecium sp. rio ete_alg 3vii", "citellus musticus", "anatidae", "escala vestjensi",
    "pogona vitticeps: feces", "monobia quadriens", "schlechtendalia chinensis", "gelatinous salpe",
    "gromphadorhina grandidieri", "folsomia candida", "dipetalonema caudispina", "crawfish",
    "scaptotrigona depilis", "dermacentor variabilis adult female tick", "wild goat",
    "myotis lucifugus", "bombus impatiens", "litomosoides sigmodontis", "callodictya krueperi",
    "macrotermes natalensis", "channa argus (cantor) snakehead fish", "new zealand white rabbit",
    "california sea lion", "ostrich", "ornithodoros hermsi (tick)", "neotermes castaneus",
    "salmo salar (atlantic salmon)", "warmwater fish", "channel catfish from pond", "riptortus pedestris",
    "bactericera cockerelli", "ochotona pallasi", "testudo hermanni hermanni", "soft tick",
    "ochotona curzoniae", "bos grunniens (yak)", "daphnia cucullata (water flea)",
    "puntigrus tetrazona", "procambarus clarkii", "flies", "procamarus clarkia", "bodo saltans",
    "draeculacephala minerva (green sharpshooter)", "larus sp. (gull)", "animals",
    "antarctic penguin", "planococcus citri (mealybug)", "buteo magnirostris", "hare",
    "amblyomma sculptum", "meimuna oshimensis", "polyplax serrata", "threskiornis spinicollis",
    "bugula neritina", "crocidura russula", "sitophilus oryzae", "epithemia turgida",
    "eulemur macaco", "cockroach", "salmon", "ctenarytaina eucalypti", "perna canaliculus",
    "wildlife", "haematopota sp. (horse fly parasite)", "oncorhynchus kisutch (coho salmon)",
    "salganea taiwanensis taiwanensis", "pygoscelis papua", "camponotus chromaiodes; colony 640",
    "macrosiphum euphorbiae", "thais luteostoma", "epinephelus coioides", "pygoscelis papua",
    "euplotes petzi", "putative tunicate", "rice leafhopper", "balta sp.", "colobopsis nipponica",
    "brugia malayi strain trs", "panesthia angustipennis", "chonosia crassipennis",
    "carbrunneria paramaxi", "oliarus filicicola (zimmerman 1945)", "spermophilus pygmaeus",
    "mimachlamys nobilis", "homalodisca vitripennis (glassy-winged sharpshooter)",
    "mactra quadrangularis", "spheniscus humboldti", "yellow catfish", "honeybee", "euplotes sp.",
    "enhydra lutris nereis", "girella melanichthys", "panthera tigris amoyensis facel",
    "passer domesticus", "rhyparobia maderae", "strigomonas culicis", "mimachlamys nobilis",
    "larus michahellis", "camponotus nipponensis", "equus ferus", "turkeys", "myotis velifer",
    "dermacentor variabilis", "camponotus (colobopsis) obliquus", "schizaphis graminum",
    "panesthia angustipennis yayeyamensis", "sparus aurata", "dysmicoccus sylvarum",
    "monobia quadridens", "mesoclemmys nasuta", "llaveia axin axin", "marine sponge",
    "spheniscus magellanicus", "vagitanus terminalis", "caprine (goat)", "amblyomma ovale",
    "macrosteles sp. nr. severini", "tiphiidae", "mandarin fish", "small ruminants",
    "black-faced spoonbill", "murine", "lateolabrax japonicus", "meimuna opalifera",
    "allacta australiensis", "sipunculus nudus l.", "aedes aegypti (rockefeller strain)",
    "tabanus abactor", "bactericera cockerelli (potato psyllid)", "dolichoderus sibiricus",
    "zootermopsis nevadensis", "platypleura kaempferi", "ciconia ciconia", "orthaga achatina",
    "salmonid fish", "mice", "clastoptera arizonana (arizona spittlebug)", "malayan pangolin",
    "calcareous sponge clathrina clathrus", "tabanus", "varanus sp.", "tabanus sp.",
    "stachyamoeba lipophora atcc 50324", "spumella elongata", "pachypsylla sp. 'celtidis'",
    "african ground squirrel", "phocoena phocoena (harbor porpoise)", "rabbit tick",
    "tursiops truncatus (bottlenose dolphin)", "bemisia tabaci q", "grasshopper",
    "aphis gossypii", "drosophila simulans riverside", "japanese macaque (macaca fuscata)",
    "aedes albopictus", "penaeus plebejus", "apis mellifera (western honey bee)",
    "rock hyrax", "chrysomya megacephala (blowfly)", "midday gerbil", "iguana",
    "leptinotarsa decemlineata (colorado potato beetle)", "caracal caracal",
    "drosophila ananassae florida (14024_0371.12)", "tilapia (oreochromis niloticus)",
    "lophura swinhoei", "marine sponge pandaros acanthifolium", "sea lion",
    "finch", "anopheles sinensis", "cicindelidae", "hirondellea gigas",
    "heteropsylla cubana", "syrphidae", "salvelinus namaycush (lean lake trout)",
    "heteropsylla texana", "snail", "wild rodents", "sponge", "hymeniacidon perleve",
    "dirofilaria (dirofilaria) immitis", "graphocephala atropunctata", "planiliza haematocheilus",
    "mogannia minuta", "microtus fortis", "chauliognathus", "rhynchophorus ferrugineus",
    "macrobrachium rosenbergii", "holotrichia parallela", "steno bredanensis",
    "bemisia tabaci", "camel", "tabanus", "finespotted flounder", "ochotona curzoniae (pika)",
    "cunner fish", "aspidiotus nerii bouche", "crocodylus niloticus", "seal", "cantharis",
    "bluespotted ribbontail ray (taeniura lymma)", "hemlock woolly adelgid", "psyllid",
    "bathymodiolus thermophilus", "macrosteles quadrilineatus", "vulpes vulpes",
    "macrosteles quadripunctulatus", "ursus americanus", "tibetan wild ass", "peacock",
    "diuraphis noxia", "koala", "black soldier fly", "gyna capucina", "guinea fowl",
    "oreochromis niloticus (tilapia)", "hyalessa maculaticollis", "kosemia yezoensis",
    "netuma thalassina", "yezoterpnosia vacua", "mirounga angustirostris (northern elephant seal)",
    "eriocheir sinensis", "calyptogena fausta", "sea squirt", "paratemnopteryx sp.",
    "passer hispaniolensis", "caenorhabditis elegans my316", "pika", "cinereous vulture",
    "pantholops hodgsonii (tibetan antelope)", "euphyllodromia sp.", "hilsa fish",
    "bemisia tabaci (biotype b)", "liposcelis bostrychophila", "microtus arvalis",
    "paulinella", "ranissus scytha", "cicadellidae", "tetraponera penzigi",
    "tettigades undata", "ctenarytaina spatulata", "tettigades limbata", "sogatella furcifera",
    "gut of protaetia brevitarsis seulensis (larva)", "entylia carinata", "ectobius sp.",
    "oncorhynchus mykiss (rainbow trout)", "drosophila mauritiana", "meriones meridianus",
    "ctenopharyngodon idella (grass carp)", "siniperca chuatsi", "oncorhynchus tshawytscha",
    "tabanus nigrovittatus (horsefly)", "reticulitermes speratus", "curculio caryae",
    "rabbit tick", "rhopalodia gibberula", "ornithodoros parkeri", "zalophus californianus",
    "auritibicen bihamatus", "prairie dog", "protaetia brevitarsis seulensis larva",
    "diestrammena coreana (camel cricket)", "sycon capricorn", "euscepes postfasciatus",
    "mediastinia sp.", "cosmozosteria sp.", "amblyomma americanum", "platyzosteria sp.",
    "cellana toreuma", "ring-tailed lemur", "testudinidae", "photinus",
    "haemaphysalis concinna", "apis cerana", "heterodera glycines", "pyxicephalus edulis",
    "trialeurodes vaporariorum", "ciconia boyciana (oriental stork)", "tettigades chilensis",
    "cryptotympana atrata", "bemisia tabaci china 1", "kentomonas sorsogonicus",
    "animal", "lizard", "wildebeest", "urostylis westwoodii", "steinernema monticolum",
    "phoxinus keumkang", "sturgeon", "heterorhabditis bacteriophora", "gerbil",
    "ixodes ovatus", "otter", "odobenus sp.", "pachyrhynchus infernalis", "panorpa helena",
    "epinephelus aeneus", "sander vitreus", "sheep (ovis)", "anser albifrons",
    "pantholops hodgsonii", "cats", "capra aegagrus hircus", "cimex lectularius jesc",
    "spermophilus dauricus", "diceroprocta semicincta (cicada)", "frog",
    "donacia thalassina", "ixodes persulcatus", "tegillarca granosa", "wild bird",
    "edessa sp.", "nematode", "seriola lalandi", "marmota himalayana", "egret",
    "tragelaphus strepsiceros", "apis mellifera", "botryllus sp.", "oncorhynchus kisutch",
    "meleagris gallopavo", "blue sheep", "shellfish", "insect", "chelymorpha alternans",
    "neophocaena asiaeorientalis", "penaeus japonicus", "alligator", "macroplea mutica",
    "tibetan antelope", "penaeus vannamei (whiteleg shrimp)", "larus sp.", "protozoa",
    "donacia provostii", "donacia marginata", "ictalurus punctatus", "amblyomma cajennense",
    "argas persicus", "mink", "cistudinella sp.", "ischnocodia annulus", "bigeye tuna",
    "drosophila melanogaster oregon-r-modencode", "sablefish", "oncorhynchus mykiss",
    "atlantic salmon", "shrimp", "amblyomma variegatum (cattle tick)", "nile tilapia",
    "tilapia fish", "migratory bird", "south china tiger", "crocodile lizard",
    "arvelius albopunctatus", "crassostrea gigas", "phoca largha", "misgurnus anguillicaudatus",
    "xenopus laevis", "trigona sp.", "buteo jamaicensis", "aphis fabae", "eclectus roratus",
    "equus kiang", "phractocephalus hemioliopterus", "black-collared starling", "marmot",
    "cryptocercus clevelandi", "wild chukar", "alvinocaris longirostris",
    "penaeus vannamei", "donacia versicolorea", "nipponaphis monzeni", "salmo salar",
    "amazona sp.", "periplaneta americana", "mizuhopecten yessoensis", "haliotis discus hannai",
    "plateau pika", "silkworm", "gentoo penguin", "discomorpha sp.", "tibetan gazelle",
    "chroicocephalus novaehollandiae (australian silver gull chick)", "trachinotus ovatus",
    "crow", "cassida sp.", "pelodiscus sinensis", "hirudo verbana", "procambus clarkii",
    "galaxea fascicularis (stony coral)", "bivalve mollusk", "ondara zibethicus",
    "aphis glycines", "mus musculus", "sebastes schlegeli", "mus musculus subsp. domesticus",
    "donacia semicuprea", "tuberolachnus salignus", "hippopotamus amphibius", "jellyfish",
    "paralichthys olivaceus", "bos mutus", "pheasant", "charidotella sexpunctata",
    "scylla serrata", "euphausia superba", "canis", "adelges kitamiensis",
    "neoaliturus tenellus", "crucian carp", "ictalurus furcatus", "mercenaria mercenaria",
    "physopelta gutta", "agroiconota sp.", "pentalonia nigronervosa", "bison bison",
    "apterostigma dentigerum", "urocitellus undulatus", "mussels", "artemisaphis artemisicola",
    "aphis nerii", "eel", "ixodes ricinus (tick)", "acromis sparsa", "carposina sasakii",
    "macroplea appendiculata", "donacia clavipes", "connochaetes taurinus", "perinereis linea",
    "trichonympha agilis", "wild boar", "pseudotrichonympha sp.", "flounder",
    "drosophila melanogaster", "dicentrarchus labrax", "rainbow trout", "bank vole (clethrionomys glareolus)",
    "donacia dentata", "oyster larvae", "myocastor coypus", "megalobrama amblycephala",
    "honey bee", "channa argus (snakehead fish)", "crocodylus siamensis", "forficula auricularia",
    "ixodes scapularis", "gut of wasp", "american cockroach", "trypoxylus dichotomus",
    "parrot", "paper wasp", "marmota", "scophthalmus maximus", "locusta migratoria",
    "oreochromis niloticus", "ornithodoros hermsi", "oecophylla smaragdina", "euschistus servus",
    "ixodes ricinus", "cassida viridis", "cassida vibex", "cynoglossus semilaevis",
    "phakellia ventilabrum", "phoca vitulina", "megacopta punctatissima", "invertebrates",
    "acyrthosiphon pisum", "donacia vulgaris", "mouse", "mya arenaria oonogai makiyama",
    "nasonia vitripennis", "cryptocercus punctulatus", "porcellio scaber", "oreochromis sp.",
    "schizaphis graminum biotype i", "oreochromis", "cage-cultured red drum", "coral",
    "ornithodoros sonrai", "fish", "amblyomma longirostre", "ceratitis capitata",
    "gorilla", "hydra vulgaris aep", "amblyomma neumanni", "plateumaris rustica", "oryctes gigas",
    "crab", "honeybees", "aratinga solstitialis", "correlophus ciliatus", "struthio camelus",
    "syngnathus typhle", "nezara viridula (cotton pathogen vector; southern green stink bug)",
    "deinagkistrodon acutus", "aphis helianthi", "anguilla anguilla", "halyomorpha halys",
    "antho dichotoma", "grass carp", "bighead carp", "trichechus manatus", "anguilla japonica",
    "pthirus gorillae", "marmota baibacina", "ruddy shelduck", "marmota sibirica",
    "protaetia brevitarsis seulensis", "cancer pagurus", "lepus brachyurus", "lepus europaeus",
    "bactrocera oleae", "rattus", "hyadaphis tataricae", "myzus persicae", "alces alces",
    "flatfish", "toothfish", "bankia setacea", "litopenaeus vannamei", "naegleria",
    "acyrthosiphon kondoi", "sebastes schlegelii", "penaeus vannamei (shrimp)", "ovis",
    "macrosiphum gaurae", "turbot", "turbo", "adelges lariciatus", "adelges cooleyi",
    "pseudotrichonympha grassii (protist) in the gut of the termite coptotermes formosanus",
    "cerambycidae sp.", "pediculus humanus corporis", "ixodes pacificus", "donacia proxima",
    "kangaroo", "cryptopsaras couesii", "anas platyrhynchos", "penaeus setiferus",
    "apis mellifera mellifera", "culicoides impunctatus", "larimichthys crocea",
    "camponotus chromaiodes", "caenorhabditis elegans", "murgantia histrionica",
    "parakeet", "geronticus eremita", "melanaphis sacchari", "dendrolimus ibiricus",
    "drosophila melanogaster oregon", "allacta bimaculata", "peafowl", "melopsittacus undulatus",
    "homalodisca coagulata (glassy-winged sharpshooter)", "pachypsylla venusta",
    "cinara pseudotaxifoliae", "cuttlefish",
    "myzus persicae (green peach aphid)", "macrotermes barneyi", "glossina morsitans morsitans",
    "nauphoeta cinerea", "delphinapterus leucas (beluga whale)", "sea urchin", "cinara piceae",
    "oryctolagus cuniculus", "corvus brachyrhynchos", "drosophila melanogaster oregon-r modencode",
    "apterostigma", "helicoverpa armigera", "seabass", "rabbit", "plateumaris braccata", "rat",
    "odocoileus virginianus", "peromyscus leucopus", "cinara cedri", "brachycaudus cardui",
    "plateumaris consimilis", "acanthamoeba polyphaga hn-3", "perameles bougainville",
    "hyalomma aegyptium", "zebra", "apsterostigma", "baizongia pistaciae", "pediculus schaeffi",
    "otospermophilus beecheyi (california ground squirrel)", "apodemus agrarius", "anguilla japonica",
    "paralichthys olivaceus (flounder)", "dermacentor andersoni", "donacia tomentosa", "lacertilia",
    "papio papio", "chroicocephalus novaehollandiae", "ochotona curzoniae (plateau pika)",
    "anomalops katoptron", "gull", "japanese eel", "salvelinus fontinalis", "blattella germanica",
    "pteropus livingstonii", "pteropus poliocephalus", "bat", "skate", "red kangaroo",
    "moschus berezovskii", "macaca silenus", "panulirus ornatus", "pavona duerdeni",
    "odontobutis platycephala", "lutjanus guttatus (rose snapper)", "kentomonas sorsogonicus",
    "ostrea edulis (flat oyster)", "beluga whale", "shinkaia crosnieri", "haematopota pluvialis",
    "pelteobagrus fulvidraco", "aquila chrysaetos", "graptopsaltria bimaculata",
    "vultur gryphus", "penaeus monodon", "tanakia koreensis", "dicentrarchus labrax",
    "zophobas atratus", "pigeon", "turtle", "penaeus (litopenaeus) vannamei (whiteleg shrimp)",
    "donacia cinerea", "fulmars", "brevicoryne brassicae", "acanthamoeba", "blaberus giganteus",
    "arctic char", "microlophium carnosum", "photinus pyralis", "ciconia boyciana",
    "macaca mulatta", "thelaxes californica", "bactrocera dorsalis", "apostichopus japonicus",
    "ovis aries", "leptinotarsa decemlineata", "phyllophaga sp.", "bombyx mori",
    "pthirus pubis", "manis javanica (pangolin)", "sipalinus gigas", "geodia barretti",
    "macrosiphoniella sanborni", "muscaphis stroyani", "aphis urticata", "cotinis nitida",
    "crassostrea gigas (pacific oyster)", "crassostrea virginica", "ellychnia corrusca",
    "glossina brevipalpis", "marine sponge lissodendoryx isodictyalis in the bahamas",
    "catfish", "morone chrysops x morone saxatilis", "pediculus humanus", "rattus norvegicus",
    "bothriocroton concolor", "hyperomyzus lactucae", "trimyema compressum", "macaca fascicularis",
    "neohaemonia nigricornis", "haemaphysalis juxtakochi", "columba livia", "diaphorina citri",
    "prawn", "deer", "nebria ingens riversi", "donacia fulgens", "sitobion avenae",
    "pseudorca crassidens", "acyrthosiphon lactucae", "aphis nasturtii", "giant panda",
    "cryptocercus kyebangensis", "donacia sparganii", "oyster", "caprine", "rattus sp.",
    "japanese rhinoceros beetle larva", "hydrophilus acuminatus", "discus", "hawk",
    "acinonyx jubatus", "chaeturichthys stigmatias", "donacia cincticornis", "anas strepera",
    "pan troglodytes", "paguma larvata", "plateumaris sericea", "diplonemea", "cybister lewisianus",
    "aphis craccivora (cowpea aphid)", "calanoid copepod", "lagenorhynchus acutus",
    "ctenocephalides felis", "parachirida sp.", "nezara viridula", "euterpnosia chibensis",
    "pyropia yezoensis conchocelis", "pan paniscus", "lipaphis pseudobrassicae",
    "donacia piscatrix", "diaphorina cf. continua", "lion-tailed macaques", "penguin",
    "aphis craccivora", "donacia crassipes", "donacia bicoloricornis", "donacia simplex",
    "draeculacephala minerva", "spodoptera frugiperda", "danaus plexippus", "danio rerio",
    "anodonta arcaeformis", "mastotermes darwiniensis", "toxic alexandrium minutum",
    "coreoleuciscus splendidus", "avian", "mus musculus c57bl/6j", "tick", "rhesus macaque",
    "platycercus elegans", "solea senegalensis", "sanzinia madagascariensis volontany",
    "plateumaris pusilla", "gilthead seabream", "ixodes pacificus (western blackleg tick)",
    "coregonus clupeaformis (lake whitefish)", "wax moth", "galleria mellonella",
    "cyclopterus lumpus", "cockatiel", "melanocetus johnsonii", "acheta domesticus",
    "blatta orientalis", "seriola dumerili", "white stork", "canis latrans", "scallop",
    "drosophila neotestacea", "shrimps", "bird of prey", "ardea cinerea", "korean rockfish",
    "c57bl/6ntac mice", "phascolarctos cinereus", "neopsylla setosa", "haliclona simulans",
    "rhopalosiphum padi", "silurus asotus", "uroleucon ambrosiae", "mustela putorius furo",
    "meimuna kuroiwae", "tettigades ulnaria", "cystophora cristata", "porites pukoensis",
    "ringed seal", "holothuroidea", "tilapia", "macropus giganteus", "haliotis gigantea",
    "ailuropoda melanoleuca", "spermophilus sp.", "squirrel", "dendrolimus sibiricus", "cormorant",
    "melanozosteria sp.", "protagonista lugubris", "therea regularis", "owl", "mule",
    "hymenochaete rubiginosa", "himantormia", "veronaeopsis simplex y34", "podila verticillata",
    "fusarium oxysporum f. sp. cucumerinum", "flammulina filiformis", "rhizoctonia solani",
    "white-rot fungus phanerochaete chrysosporium", "mushroom", "flammulina velutipes",
    "cinara kochiana kochiana", "cinara curvipes", "cinara splendens",
    "dicentrarchus labrax (sea bass)", "cinara laricifoliae", "cinara strobi", "cinara pseudotsugae",
    "entomortierella parvispora", "alexandrium minutum", "cinara cf. splendens/pseudotsugae 3390",
    "mortierella parvispora", "shiraia bambusicola", "stereocaulon"]


animal_isolation_sources = [
    "myodes rufocanus", "terpios hoshinota", "fish (cod)", "chimpanzee", "a dead ark clam",
    "sea squirt", "archimandrita tessellata", "gut of tenebrio molitor", "saimiri sciureus", "silkworm feces",
    "cacatua galerita", "panthera tigris amoyensis", "geotrupidae", "paguma larvata",
    "animal (cat)", "chlorocebus sabaeus", "ailuropoda melanoleuca", "bamboo coral",
    "re-isolation from dsm 1974", "spiral shell", "miss", "mice feces", "seastar",
    "rabbit tissue (derivative of egd strain)", "bacteriocytes of mediterranean bemisia tabaci",
    "calcareous sponge clathrina clathrus", "oncorhynchus mykiss", "gelatinous salpe",
    "isolated from adult of haemaphysalis hystricis", "monobia quadriens", "pogona vitticeps: feces",
    "sturnus nigricollis", "rabbit tick", "spoon worm", "meleagris", "sparus aurata", "ornithodoros turicata",
    "whole fish sold as black tilapia; muscle", "shellfish: shrimp", "wild bird feces",
    "onchrhynchus mykiss", "farmed eel", "chitin", "environmental-quail", "animal tissue",
    "panthera tigris amoyensis facel", "coleoptera cantharidae", "clams", "silurus asotus",
    "marmot respiratory tract", "cockatoo feces", "bee gut", "animal feces", "cuttlefish",
    "myotis velifer", "gentoo penguin", "insect gut", "female abdomens", "oysters", "nigg3",
    "grey headed albatross", "shark", "antarctic penguin faecal samples", "hindgut of honeybee",
    "odontobutis interrupta", "whole live fish sold as snakehead-haruan; muscle", "dsm1974",
    "isolated from bird (puffinus tenuirostris)", "atta cephalotes fungus garden", "dalbulus maidis",
    "giant panda feces", "gut of larva", "gill tissue", "animal - avian-environment swab",
    "plume", "sea-bass", "faeces of the tibetan antelope", "mosquito", "shell egg",
    "anodonta arcaeformis", "ephydatia sp.", "sponge", "unspecified insect", "small ruminants",
    "whiteflies bacteriome", "feces of tibetan antelope", "reserve speciale ambohitantely",
    "tissue from oliver flounder", "oyster/environment", "wild yak dung",
    "coelomic fluid of a sand dollar", "animal - avian", "mouse", "midday gerbil", "gut of limpet",
    "drosophila melanogaster", "isolated from larva of dermacentor taiwanensis", "epinephelus coioiaes",
    "rodent", "oyster", "moribund muskrat (ondatra zibethica)", "muskrat spleen", "flamingo",
    "exterior surface of the shell of an abalone sold in a fish market in tokyo; japan",
    "isolated from nymph of haemaphysalis hystricis", "tunicate", "termite gut", "thais luteostoma",
    "murine model", "mus musculus", "liposcelis bostrychophila", "hilsa fish", "hill myna",
    "isolated from haemaphysalis hystricis", "squid", "abalone haliotis discus hannai",
    "feces of an oriental stork (ciconia boyciana)", "reptile", "tortoise tubercle",
    "crayfish", "iguana iguana", "crustacean shell", "malayan pangolin", "yak feces",
    "animal", "wild rat", "coho salmon", "fruit fly", "diaphorina citri", "south china tiger",
    "iguana", "bat", "wild yak feces", "white-crowned sparrow", "african elephants",
    "efb-infected honey bee colony", "turban shell", "mussels", "rainbow trout",
    "foulbrood of honeybees", "turbot", "red fox", "panda", "neopsylla setosa",
    "marine sponge", "chamaeleonidae", "ixodes spinipalpis", "barnacle at wood pile-on",
    "atlantic salmon", "manis javanica", "a feces sample of migratory birds origin",
    "tilapia", "wild pig; fecal", "shellfish", "diarrheal snake diarrheal snake in hunan",
    "diseased shrimp", "wasp honeycombs", "eclectus roratus feces", "frog",
    "insect larvae", "tissue; animal", "yellowtail", "snake", "hydrophilus acuminatus",
    "crassostrea gigas", "pigeon", "gill tissue of bathymodiolus japonicus", "mice (c57bl/6)",
    "meleagris gallopavo", "sea bass", "mouse gut", "avian", "sockeye salmon", "scallop larva",
    "diseased labeo rohita fish", "mealworm", "feces; gull (larus spp.)", "nematode community",
    "musca domestica", "atlantic cod", "siniperca scherzeri", "red-breasted parakeet",
    "midgut crypts of stink bug togo hemipterus", "konosirus punctatus", "heterodon nasicus",
    "crushed cell of callyspongia sp.", "tuna", "housefly", "bison", "coral", "catfish",
    "mactra veneriformis", "black-collared starling", "white-lipped deer", "elk droppings",
    "shrimp", "black-headed gull", "psittacus erithacus feces", "salmon", "silkworm excrement",
    "separated from the corpses of silkworms that had died due to bb natural infection in daiyue district; taian city; shandong province; china.", "rabbit feces", "jellyfish", "animal faecal matter",
    "rabbit", "animal hide", "tuna scrape; yellowfin", "myocastor coypus", "feces of bats",
    "cybister lewisianus", "cybister brevis", "intestinal contents of plateau pika", "vole",
    "brook charr", "penaeus japonicus",
    "aphis craccivora (cowpea aphid) on robinia pseudoacacia (locust)", "insects",
    "swiss alpine ibex feces", "psittacus erithacus", "penaeus vannamei (whiteleg shrimp)",
    "bearded dragon", "wild bird", "rattus rattus", "hepatopancreas", "pagrus major",
    "neophocaena phocaenoides", "loach", "unidentified actinians", "fish", "animal organ",
    "waxworms gut", "flea", "intestinal contents of termite nasutitermes nigriceps",
    "muskrat during outbreak", "invertebrate gut", "cockless", "herring gull; cloacal swab",
    "gromphadorhina portentosa cockroaches", "mushroom substrate", "brewery yeast",
    "fungal mycelia of mortierella elongata fmr23-6",
    "white-rot fungus phanerochaete", "hyphae", "white-rot fungus phanerochaete chrysosporium",
    "fungal hypha of mortierella elongata fmr23-6", "entomortierella parvispora"]


plant_hosts = [
    "gynura procumbens", "iris germanica", "paulownia", "chrysanthemum x morifolium",
    "mallotus japonicus", "alnus rubra", "fagus sylvatica", "alliaria petiolata", "glycyrrhiza uralensis",
    "ziziphus jujuba mill.", "abies koreana", "populus tomentosa", "vaccinium",
    "valeriana jatamansi jones", "tithonia diversifolia", "lycium barbarum l.", "lotus japonicus",
    "phalaris paradoxa", "carnation", "calystegia soldanella", "picochloruma", "xanthium sibiricum",
    "tripterygium wilfordii", "corylus avellana", "oryza longistaminata", "plants",
    "phyllosphere of grasses", "melaleuca quinquenervia", "trachydiscus minutus",
    "brasillia sp.", "cinnamomum camphora", "eutrema wasabi", "parthenium hysterophorus",
    "phaseolus microcarpus", "stevia rebaudiana bertoni", "new guinea impatiens",
    "populus euphratica oliv", "bruguiera gymnorhiza", "polypogon monspeliensis",
    "carex pumila", "adenophora trachelioides maxim", "leaves", "hypericum perforatum",
    "olive knot (olea europaea)", "basil", "astragalus pelecinus l.", "rosa rugosa",
    "linum austriacum ssp. austriacum", "leersia hexandra", "pasture gramineae",
    "pistacia chinensis bunge", "tamarix", "rehmannia glutinosa", "luffa", "poplar",
    "sagittaria sagittifolia", "acacia acuminata", "malva verticillata", "euterpe oleracea",
    "oryza latifolia", "black mangrove", "jatropha curcas l.", "bidens sp.", "pine wood",
    "cenchrus macrourus", "botryococcus braunii", "oak tree", "dactylis glomerata",
    "ulva fenestrata (green alga)", "calystegia hederacea", "dioscorea sansibarensis",
    "quercus sp.", "euonymus japonicus", "quercus", "crassulaceae", "tectona grandis",
    "tripsacum laxum", "paullinia cupana", "phyllanthus urinaria l", "fetuca",
    "adenophora trachelioides maxim.",
    "suaeda salsa", "cladonia borealis", "medicago truncatula", "ulva pertusa (algae)",
    "allamanda cathartica", "astragalus pelecinus", "lebeckia ambigua", "lotus corniculatus",
    "trifolium repens", "curcuma aromatica", "curcuma wenyujin y.h. chen et c. ling",
    "musa balbisiana cultivar kepok", "lemna trisulca", "pelargonium capitatum",
    "datisca glomerata", "beta vulgaris", "robinia pseudoacacia", "ficus religiosa l.",
    "sporobolus anglicus", "casuarina equisetifolia", "mulberry", "ficus benjamina",
    "leontopodium alpinum", "tanacetum vulgare", "catharanthus roseus", "stereocaulon sp.",
    "oxytropis pumilio", "malus sylvestris", "angelica sinensis dlies", "mountain ginseng",
    "agave americana l.", "stachytarpheta glabra", "elymus tsukushiensis", "eucalyptus",
    "androsace koso-poljanskii", "bromus inermis", "peperomia dindygulensis miq.",
    "trifolium pratense", "lemna minor", "mimosa affinis", "salsola stocksii", "mimosa",
    "biserrula pelecinus l.", "rhizoma kaempferiae", "masson pine", "populus x jackii",
    "nerium oleander", "biserrula pelecinus", "lichen", "miscanthus giganteus", "mimosa scabrella",
    "rosa sp.", "malus sieversii", "onobrychis viciifolia", "carex sp. (sedge blades)",
    "flower", "mimosa flocculosa", "trifolium uniflorum", "eucalypti of eucalyptus", "wild cotton",
    "urochloa reptans", "cotinus coggygria", "kandelia candel (mangrove)", "mangrove",
    "trifolium spumosum l. (annual mediterranean clovers)", "cardamine cordifolia",
    "medicago orbicularis", "arabidopsis", "vachellia farnesiana", "arabidopsis thaliana",
    "physcomitrium patens",
    "physcomitrella patens", "cryptomeria japonica var. sinensis", "festuca arundinacea", "oxytropis triphylla",
    "usnea", "rubus sp.", "melilotus officinalis", "himalayan blackberry", "wild rice species",
    "achillea ptarmica", "muscari", "lespedeza cuneata", "nerium oleander", "algae",
    "tanacetum vulgare", "eucommia ulmoides", "clematis", "quercus rubra", "dendrobium officinale",
    "leiosporoceros dussii", "oxytropis triphylla", "vavilovia formosa", "oxytropis kamtschatica",
    "zantedeschia aethiopica", "rosa", "prunus cerasifera", "medicago", "phalaenopsis sp. (orchid)",
    "leontopodium nivale", "nusuttodinium aeruginosum", "prosopis cineraria", "clivia",
    "oryza glumipatula", "acacia farnesiana", "catalpa", "populus", "jatropha curcas",
    "populus alba x (p. davidiana + p. simonii) x p. tomentosa", "antarctic ice algae",
    "salvia splendens", "panax ginseng", "vicia alpestris", "phaeoceros", "trifolium sp.",
    "commiphora wightii", "alhagi sparsifolia shap.", "acacia farnesiana", "hyacinthus sp.",
    "lobaria pulmonaria thallus", "lotus", "erigeron annuus l. pers", "sida hermaphrodita",
    "hyacinthus orientalis", "sporobolus anglicus", "malus prunifolia (crab apple)",
    "tephrosia apollinea", "vitis", "echinacea purpurea", "phyla canescens", "ripariosida hermaphrodita",
    "aeschynomene fluminensis", "vavilovia formosa", "prunus sp.", "plant tissues",
    "oenothera speciosa", "acaciella angustissima", "brassica rapa subsp. chinensis",
    "catharanthus roseus", "dracaena sanderiana", "parthenium argentatum gray (guayule shrubs)",
    "medicago arborea", "medicago lupulina l.", "dongxiang wild rice", "tephrosia purpurea subsp. apollinea",
    "euonymus sp.", "eucalyptus grandis", "broussonetia papyrifera", "lichen stereocaulon sp.",
    "lotus sp.", "weed", "thlaspi arvense", "grass", "quercus castaneifolia", "miscanthus x giganteus",
    "blasia pusilla", "calliandra grandiflora", "brittle root", "digitaria eriantha",
    "vavilovia formosa", "bossiaea ensata", "liriodendron tulipifera", "pine tree",
    "plant", "sesame", "salix sp.", "hibiscus", "hibiscus rosa-sinensis"]



plant_isolation_sources = [
    "mulberry mosaic dwarf leaf", "surface-sterilized vetiver roots", "arctic grass",
    "white-flowered calla lily", "luffa aegyptiaca", "root rhizomes", "tree", "sedge blades",
    "kumarahou flower", "the root nodule of glycyrrhiza uralensis", "saccharum officinarum",
    "flower of dendranthema zawadskii", "avicennia marina", "moss sample of penguin habitat",
    "bruguiera gymnorhiza", "calla lily", "mountain ginseng", "phyllosphere of grasses",
    "4-year-old roots of korean ginseng", "new guinea impatiens", "mangrove", "algal mat",
    "sida hermaphrodita", "roots of desert plants", "pphyllosphere of grasses",
    "arabidopsis thaliana seedling from surface sterilized seed", "plant matter", "grass",
    "flower of forsythia koreana", "algal phycosphere", "lotus corniculatus nodule",
    "malus sylvestris", "morning glory", "dongxiang wild rice", "napier grass", "carex pumila",
    "flower of rhododendron schlippenbachii", "root of codonopsis pilosula", "root nodules",
    "chrysochromulina tobin phycosphere", "algae", "artificially infected catharanthus",
    "leaf spots", "poplar tree gall", "gardenia fruit", "plant root", "bark of cypress", "fig leaf",
    "peltigera membranacea thallus", "angelica gigas nakai root surface", "plant root nodule",
    "galega orientalis root nodule", "seaweed", "xylem", "red seaweed", "macrocystis pyrifera",
    "plant", "broussonetia papyrifera", "lichen", "lowbush blueberries", "root", "roots",
    "seaweed; enteromorpha linza (coastal marine sediments); aburatsubo inlet; australia",
    "flower", "leaf of fig tree", "grass from the park in a vetenary clinic", "halobiotic reed",
    "plant roots", "galega officinalis root nodule", "lichen thallus", "ginseng", "gladiolus",
    "sesbania spp. root nodule", "plant leaf", "origanum marjorana", "malus sp.", "algae",
    "flower of shittah tree", "himalayan blackberry", "endosphere environment",
    "inner tissues of halophyte limonium sinense (girard) kuntze", "leaf", "hiacynth plant",
    "root nodule of sesbania cannabina", "wood", "necrotic oak lesion", "tree root",
    "fallen leaves of virgin forest", "flower of chiness redbud", "calliandra haematocephalus",
    "flower of rhododendron sclippenbachii", "plant xylem"]

terrestrial_isolation_sources = [
    "top of mount daeam; south korea", "yuncheng salt lake", "stratosphere", "envo:01001711",
    "nature biofilm", "red biofilm", "molybdenum mine", "biofilms", "natural biofilms", "atcc 43099",
    "pieces of wood at the bottom of a cave", "environment monitoring (air)", "natural biofilm",
    "2-km-deep aquifer", "aquifer", "richmond mine", "air sample", "debris flow", "biofilm material",
    "crystal cave (limestone)", "karstic limestone", "granitic rock aquifer at 600 m depth",
    "terrestrial subsurface brine", "a deep subsurface coal seam", "bentonite", "mineral material",
    "obsidian pool; yellowstone national park", "water from a salt pit", "mat", "taupo volcanic zone",
    "stomatolites grown in far red light (720nm)", "surface of island", "surface brine",
    "microbial mat; hamelin pool", "tar pits", "salt field", "salt particles", "limestone cave wall",
    "sulfur enrichment", "environmental sample", "outdoor air", "sulfide deposits",
    "salt crystallizer of the great rann of kutch", "alkaline pod", "gomso solar saltern",
    "isolated from slime streamers and attached to pyrite surfaces at a sulfide ore body",
    "2.8-km deep subsurface aquifer", "subsurface mine microbial mat", "hot spings runoff",
    "64.1 degree centigrade; ph 3.72 in naghaso; the philippines", "cave entrance",
    "64.1 degree centigrade; ph 3.72", "rock from karst cave", "dead sea", "shale",
    "solfataric thermal field close to moutnovsky volcano", "halite rock", "qaidam basin",
    "iron hydroxide deposits", "solar salt", "air", "elton hypersaline lake",
    "brazilian saline-alkaline lake", "desert", "saltpan", "hyperthermophilic compost",
    "coal bed", "moderate hot spring", "salt marsh", "alkaline pool submerged anode electrode",
    "kulunda steppe hypersaline lake", "hypersaline environment", "saline lake",
    "chronically low temperature and dry polar region", "thermophilic environment",
    "solar saltern", "dust", "beach sand", "salar de atacama; atacama desert", "beach",
    "weathered rock sample", "antimony mine", "sand", "weathered tuff", "salt crystallizer",
    "moist arsenopyrite (feass)-containing rock taken from a mine tunnel approximately 300 m below the ground in the granites gold mine", "microbial mat", "subsurface rock", "coal", "microbial mat material",
    "salt crystallizer of little rann of kutch", "stalactite biofilm", "saturated brine",
    "inside the caves of drach", "non-purified solar salt", "gold-copper mine",
    "antimony vein of nakase mine",
    "xinjiang aibi salt lake", "carbonated precipitates", "soda-saline lake", "microbial mat/biofilm",
    "shar-burdin hypersaline soda lake", "marine solar saltern brine", "soda lake magadi",
    "rock (envo:00001995)", "geothermal reservoir", "rocky sand", "solar saltern in gomso bay",
    "lava", "brine", "saltern", "baengnokdam summit crater area; mt. halla", "saltern pond",
    "taibei marine solar saltern near lianyungang city", "orthoquartzite cave surface", "cave",
    "the bange salt-alkaline lake in tibet", "the surfaces of weathered potassic trachyte",
    "deep subseafloor coal bed", "deep subsurface anoxic brine", "chaqia salt lake",
    "phototrophic microbial mat in hot lake; a shallow mgso4 dominated salt lake",
    "geothermal isolate", "salt crust", "coal spoil heap", "salt lake", "acidic salty water",
    "flat; laminated microbial mat in a salt marsh", "solar saltern of 19% salinity",
    "rock (envo:00001995)", "saline saltern", "hypersaline lake", "yates shaft; surf",
    "thiodendron' bacterial sulfur mat from mineral sulfide spring"]

soil_hosts = ["decaying wood"]

soil_isolation_sources = [
    "trebel valley fen", "rhizobacterium", "heilongjiang province", "solo", ## not an error: portuguese for soil.
    "rice paddy field", "mixed forest with bamboo", "mixed forest near oak tree",
    "envo:00001998",
    "decomposing algal scum", "rice paddy", "envo:00005769", "pasture", "peanut land",
    "decaying vegetation", "corn stalk residue compost", "environment - litter/manure",
    "decaying manure", "tomato rhisosphere", "antonio's farm; antonio rd.", "rhizophere",
    "agricultural site", "equine field isolate", "kaolin clays", "banks of congaree river",
    "tomato field", "agricultural fields", "cotton field", "wheat field", "pasture", "surface compost",
    "isolated from the rice field", "sugar beet rhizoshpere", "mango ochard",
    "rice fields", "sugarcane field", "less and more altered tuff", "organic material",
    "agricultural field", "tundra wetland", "permafrost", "compost", "mixed sand sample",
    "sphagnum peat from the bog obukhovskoe (acidic wetland)", "norway spruce forest humus",
    "paddy field; sungai manik; malaysia", "dune grassland", "decaying wood", "forest soi",
    "rice fields", "sphagnum peat", "soil around hot spring", "dirt", "ginseng field",
    "sphagnum bog", "fertilizer", "ancient permafrost from mammoth", "montane grasslands",
    "solar salt farm", "enriched culture of compost", "farmland", "permafrost; kolyma lowland",
    "siberian permafrost", "particulate matter", "pine forest", "mushroom compost",
    "manure compost", "angelo meadow plot 1; 20cm depth; 2 days after second rain event (91mm)",
    "composted cattle manure", "agricultural waste material", "wood decay material",
    "long-term organic manure fertilized", "compost sample from farm", "rice field",
    "Microscale soil grain", "ermafrost region of qilian mountains"]

agri_hosts = [
    "dioscorea alata", "cichorium endivia", "oryza sp.", "nicotiana tabacum", "potato_root",
    "brassica napus subsp. oleifera", "tobacco_roots_ghana", "dragon fruit tree", "rape", "chili",
    "vitis labrusca x vitis vinifera cultivar black queen", "coriandrum sativum", "eutrema japonicum",
    "brassica rapa subsp. oleifera", "wheat root", "pistacia vera l.", "brassica rapa",
    "ribes nigrum", "helianthus annuus", "fruit", "potato stem", "brassica oleracea var. capitata",
    "vasconcellea heilbornii", "triticum aestivum l. cv. tabasi", "surface-sterilized wheat roots",
    "cabbage", "phaseolus lunatus", "tapioca", "sweet potato", "pyrus communis",
    "coffea", "grape", "citrus spp.", "prunus avium", "capsicum sp.", "phaseolus vulgaris",
    "sugarcane", "lolium perenne", "zingiber officinale", "raspberry", "rapeseed plant",
    "phaseolus", "punica granatum", "pyrus pyrifolia var. culta", "alfalfa", "cucumis sativus",
    "cucurbita maxima", "lentil", "lathyrus sativus", "coffea arabica", "calimyrna fig",
    "solanum melongena", "vitis vinifera", "vitis vinifera l. cv. seto giants", "grapevine",
    "sakura tree", "rice", "japanese radish", "black pepper plant", "radish",
    "chinese cabbage", "kiwifruit", "pisum sativum l. (pea)", "common bean", "melon", "pepper",
    "glycine max", "glycine", "g.max", "manihot esculenta", "capsicum annuum", "cashew",
    "actinidia deliciosa", "sugarcane", "actinidia chinensis", "lettuce", "olea europaea",
    "grapefruit", "brassica juncea var. foliosa", "apple", "strawberry", "vasconcellea x heilbornii",
    "brassica oleracea var. botrytis", "camellia sinensis", "juglans regia", "coffee",
    "brassica rapa subsp. pekinensis", "humulus lupulus", "morus alba", "organic baby spinach",
    "triticum aestivum", "areca catechu", "date palm", "zea mays l.", "ginger", "plum", "cowpea",
    "potato", "wheat", "malus domestica", "brassica oleracea", "rice plant", "morus alba l.",
    "saccharum officinarum", "cotton", "soybean", "mandarin orange", "fava bean",
    "vitis vinifera l. cv. aurora black", "ziziphus mauritiana lam", "pyrus pyrifolia",
    "sesame seedling", "tobacco", "peanut", "carrot", "apple tree", "banana",
    "soybean (glycine max (l.) merrill)", "ipomoea aquatica", "elaeis guineensis",
    "panicum miliaceum", "solanum lycopersicum", "plantain", "pear", "winter wheat",
    "pyrus communis 'williams'", "fragariae ananassa", "hordeum vulgare", "persea americana",
    "actinidia", "zea mays", "turfgrass", "prunus dulcis", "glycine soja", "potato",
    "pepper plant", "sweet orange", "mangifera indica", "maize", "brassica oleracea var. capitata",
    "cicer arietinum", "pisum sativum", "actinidia deliciosa 'hayward'", "arachis hypogaea",
    "forage rape", "pear", "vigna radiata", "corn", "triticum aestivum (aestivum group)",
    "prunus cerasus (sour cherry)", "prunus dulcis", "sorghum bicolor", "phaseolus sp.",
    "allium cepa", "brassica oleracea var. capitata", "gossypium", "brassica juncea",
    "vitis vinifera cv. 'izsaki sarfeher'", "brassica rapa var. laciniifolia subvar. oblanceolata",
    "allium cepa l.", "blueberry", "pogostemon cablin", "glycine max cv. ac orford",
    "gossypium hirsutum", "solanum lycopersicoides", "vicia faba", "eggplant", "pear tree",
    "glycine max cv. jinju1", "lycopersicon esculentum", "camellia oleifera",
    "allium cepa (onion)", "brassica rapa ssp. pekinensis (chinese cabbage)",
    "musa sp.", "apium graveolens", "cucumber", "papaya", "brassica napus", "plant:cultivated mushroom",
    "triticum aestivum l.", "coffee plant", "pineapple", "citrus sinensis", "gossypium sp.",
    "musa spp."]

agri_isolation_sources = [
    "mango leaves", "tomato 'raisa f1'", "squash", "fresh alfalfa sprout", "pepper field", "alfalfa sprout",
    "isolated from stem bases of diseased rice plants", "dried rice straw", "carrot leaf", "cotton",
    "phaseolus vulgaris", "rhizoplane rice", "taro", "hazelnut", "wheatstraw", "raw almonds", "raw almond",
    "brassica oleracea var. capitata", "wheat anther", "grapevine phylloplane", "tobacco",
    "cucumber seedling substrate", "peach", "raw pecans", "cotton plant surface sterilized stem",
    "diseased cherry tree", "potato with the soft rot symptoms", "apple blossom", "rice stack",
    "sunflower seed hulls", "almond orchard", "oryza sativa root", "alfalfa silage", "raw pistachio",
    "parched beans", "cotton plant", "papaya", "rice straw", "potato plant", "silage", "organic baby spinach",
    "root endosphere of glycine max", "root surface of cotton; charles howell", "soybean seeds",
    "fruit", "canola root tip", "pepper plant", "corn silage", "lemon samples", "fire pearl white nectarine",
    "vegetable(malt)", "plant (maize/sorghum/rice)", "walnut", "hay", "rapeseed", "rice seeds",
    "the tomato root in nanchang", "inner tissues of cotton plant (gossypium sp.)", "raw peanut",
    "raw almond kernel nonpareil", "almond kernel (raw; variety nonpariel)", "raw pistachios",
    "leaves of field beans", "grape plant", "ginger", "baby spinach", "grape", "rotten onion",
    "styrian pumpkin anthrosphere", "citrus paradisi", "soy bean root", "cotton root", "raw pecan",
    "greenhouse grown plants with black leg infection", "marjoram", "affected cabbage tissue",
    "ucb-1 pistachio rootstock", "pineapple", "plants with moko disease", "mango", "chives",
    "raw almonds", "chilli", "eggplant", "cucumber", "rice leaves", "campbell early grape",
    "isolated from cucumber", "cabbage", "almond drupe", "mexican lime leaf", "apple twig",
    "phaseolus vulgaris root nodule", "cucumber 'kurazh f1'", "barley grains", "Pistacia vera l.",
    "solanum tuberosum", "peaches", "cornstalks and leaves", "maize seeds", "quinoa roots",
    "tomato roots", "pistachios; raw", "infected wheat", "sesame leaf", "spinach", "leaf tissues of rice in korea",
    "streptomyces spp. which was isolated from pyeonchang", "wheat roots", "corn root",
    "healthy tomato plant", "almond kernel (raw; variety carmel)", "the root of rice", "soybeans",
    "wheat germ", "raw peanuts", "rot potato tubers", "litchi pericarp", "soybean root nodules",
    "sugarcane root", "originally isolated from olive trees", "rice seed", "tea rhizoplane",
    "isolated from the tobacco substrate", "endorhiza sugar beet", "orange tree", "rice shoot",
    "hydroponic pots with potatoes", "ragi", "canola roots", "leaf from rome apple cultivar",
    "orange", "maize stem", "strawberry leaf tissue", "rye silage", "fodder", "stable grass silage",
    "beans", "carrot", "potato", "bean blight", "wilting pepper (c. annuum) stems in sanya",
    "corn", "japanese pear", "japanese peer", "strawberry leaf tissue", "wheat root", "ogi (red sorghum)",
    "oryza sativa", "nodules from common bean", "pisum sativum root-nodule", "grass silage",
    "noduls from roots of medicago sativa", "cilantro", "rice root", "rice leaf", "maize leaf",
    "naturally-infected soybean tissue", "soybean", "citrus orchard", "tomato", 
    "naturally-infected soybean leaf tissue", "apple tree", "banana", "cabbage seeds",
    "white nectarines fire pearl variety", "soybean nodules", "field isolate",
    "phyllosphere of musa spp. aaa cv. grand naine. leaf number 10"]

marine_hosts = ["argopecten purpuratus", "seawater", "pyropia tenera", "red alga", "jania sp. (red alga)",
                "saccharina japonica", "ecklonia cava", "microalgae", "sea", "thalassiosira profunda",
                "brown alga", "cyclotella cryptica", "posidonia oceanica (seagrass)", "marine",
                "(seagrass) posidonia oceanica",
                "trichodesmium erythraeum ims101", "asterionellopsis glacialis strain a3", "neopyropia tenera",
                "hard coral; acropora nasuta", "montipora capitata", "bathymodiolus septemdierum",
                "rhodosorus marinus", "red alga", "pyropia", "galaxea fascicularis", "hymeniacidon perlevis",
                "gracilariopsis lemaneiformis", "rotten red algae", "marine water", "tichocarpus crinitus (red alga)",
                "ulva prolifera", "macroalgae", "asparagopsis taxiformis", "stypopodium (algae)",
                "heterostera chilensis", "grateloupia sp.", "brown algae", "comscasia"]

marine_isolation_sources = [
    "shallow sea; symbosis with harmful algal bloom algae",
    "sea water (sand filtered); china; hong kong", "sea water [envo:00002149]",
    "the surface of red algae", "isolated from fouling in the littoral zone of the white sea",
    "turbid saltwater", "rotten brown algae", "yellow-gray fringes", "sea water of daya bay",
    "sea surface water", "coastal sea-ice sample at a depth of 0.8 m", "littoral zone",
    "the kelp laminaria japonica", "deep subsurface location at the piceance basin",
    "intertidal zone", "chimney sulfide", "deep-sea hydrothermal sulfide chimney",
    "brine water; atlantis ii deep; red sea; 2000 m depth", "brown alga", "deep-sea vent",
    "coastal surface water from the northern adriatic sea", "2216e", "micro-algae",
    "serpentinized peridotite in indian ocean", "sagami bay station p", "seashore", "seagrass",
    "leaf surface of zostera marina", "shallow subarine vent system at the kolbeinskey ridge",
    "3 m depth in a coral reef", "coastal lagoon", "japan trench; off sanriku",
    "sea water aquarium outflow; la jolla; california; u.s.a.", "chimney structure",
    "isolated by cultivation in sme medium under anaerobic conditions at 100a degrees c from a black smoker chimney fragment",
    "sme medium under anaerobic conditions at 100a degrees c from a black smoker chimney fragment",
    "collected at a depth of 1 meter; isolated by plating on solid media", "red alga",
    "surface water (5 m depth) of the central baltic sea", "chukchi sea", "east pacific ridge",
    "2600 m depth", "seaweed/starfish/sea water", "sea water of a tidal flat",
    "east pacific ocean ridge; black smoker", "chesapeake bay", "coral reef",
    "flow cytometry enriched sample from station aloha", "open ocean", "sulfide chimney",
    "edge of california current after nitrate enrichment and low light incubation",
    "black smoker wall; 3500 m depth", "surface of plastic from ocean", "sea shore water",
    "hydrothermal vent area derived", "deep sea", "ocean water", "sea-ice", 
    "water from the baltic sea", "envo:00000227", "sea", "envo:00002149", "deep-sea",
    "sea water", "deep-sea hydrothermal vent", "sea ice", "oceanic water",
    "cold seep", "envo:01000301", "tidal pool", "coral; primnoid", "east china sea",
    "ocean", "isolated from the deepest ocean", "surface of fucus serratus",
    "chrysochromulina tobin phycosphere", "gulf of finland", "cyanobacterial bloom",
    "glass slide place on reef flat in natural environment", "redoxcline baltic sea",
    "surface water of the southern north sea", "antarctic iceberg", "redoxcline black sea",
    "hydrothermal vent", "marker113 at caldera of axial seamount",
    "shallow tropical waters; normally from coral reef substrate", "open ocean water",
    "sulfidic waters (60 m) from the peruvian upwelling region", "shellfish hatchery",
    "bay of bengal", "pacific plankton", "sargasso sea", "brackish pond",
    "culture of roseovarius sp. tm1035", "coastal", "runway 10 reef (10-12 m)",
    "german wadden sea", "surface of pikea pinnata", "myojin knoll; izu-bonin arc",
    "microbial mat material from brackish estuary", "hydrothermal precipitates",
    "shallow sea; symbiosis with algae", "surface of macroalgae",
    "beach sand", "estuarine water [envo:01000301]", "150 km offshore",
    "cyanobacterial aggregates", "synechocystis sp. gt-l",
    "seawater", "shallow-sea hydrothermal system"]

freshwater_hosts = [
    "water (environment)", "water", "mesostigma viride", "groundwater",
    "microcystis culture", "microcystis aeruginosa", "cyclostephanos tholiformis",
    "skeletonema potamos", "cryoconite"]

freshwater_isolation_sources = [
    "usa: oklahoma; sallisaw; sallisaw creek public use area; robert s. kerr reservoir",
    "zhenhai reservoir in guangdong province", "seine river", "manawatu river", "shallow pond",
    "hyporheic zone", "oligotrophic pond", "aqueous", "reddish brown snow from a moor",
    "red snow", "geomnyoung pond", "hypertrophic pond", "drainage", "drainage ditch", "river water",
    "creek", "pond", "tianshan glacier", "environmental (pond)", "oligrotrophic pond", "stream water",
    "red cedar river in okemos; michigan; usa", "canal water", "freshwater from river",
    "glomma river", "snow", "microcystis culture", "glacial stream", "glacier", "rio negro river",
    "envo:00002011", "stream biofilm", "freshwater", "glacier ice", "glacial ice; 3519 m depth",
    "euglena gracilis from city ponds", "glacier from lahual spiti valley", "river", "groundwater",
    "pool at botanical garden; havana; cuba", "daechung reservoir", "ayakekumsalt lake", "envo:00002006"]

food_hosts = ["kimchi cabbage", "mung bean", "kombucha scoby", "seujeot", "pickle", "salted seafood",
              "pickled cabbage", "kimchi", "milk", "carrageenan", "frozen raw chicken", "garlic; cabbage",
              "beef", "koumiss", "nostoc flagelliforme", "sichuan pickle"]

food_isolation_sources = [
    "curd", "pumiao rice noodle", "sake", "lobster tail", "whey culture", "garam masala hot spice mixture",
    "brewing environment of maotai-flavor liquor", "pixian county bean paste", "pistachio",
    "kombucha green tea pellide", "pickling sauce", "radicchio", "daqu", "traditional youghurt",
    "sugar beet juice", "shirasu", "dry roasted pistachio", "dry sausage (skilandis)", "bottled beer",
    "refrigerated litopenaeus vannamei", "ganjang(korea soy sauce)", "sponge cake", "kem complex",
    "san francisco sourdough", "fermeted vegetable", "sausage", "myeolchi jeot", "daemi-jeot",
    "rye sourdough", "sake mash", "traditional greek wheat sourdough", "bulgogi", "rubing",
    "saeng-gimol meju", "rubing cheese", "instant soup", "soybean paste", "egusi seeds",
    "futsai", "koumiss", "soybean paste (chonggugjang)", "mashed potatoes", "fermented food",
    "kombucha scoby", "kefir", "cheese product", "instant pork", "pork &cabbage dumplings",
    "isolated from commercially available caspian sea yogurt", "broth", "mash", "kefir grains",
    "makgeolli (korean traditional alcoholic beverage)", "cream products", "bean sprouts",
    "light wheat beer", "sprouts", "3-day-old traditional semihard cheese", "pork", "frozen rock lobster tail",
    "moto starter of sake", "minced pork", "shrimp paste", "tahini (ground sesame seeds)",
    "isolated from salted brown alga laminaria", "shima dofu (okinawan-style tofu)",
    "salted brown alga laminaria", "stinky xiancaigeng", "sake(hatsuzoe)", "pla-chom",
    "blue cheese in wax", "tibicos symbiotic community", "pickled cabbage", "vegetative broth",
    "ready to eat mixed salad leaves (obtained from discount store)", "coleslaw", "cereal",
    "degraded sugar thick juice", "scallop", "chinese sausages", "buttermilk", "raw almond kernel butter",
    "ground turkey", "gimbap", "product-eggs-raw-whites", "traditional dairy product",
    "frozen peas", "palm wine", "natural whey culture from gruyere cheese",
    "homemade koumiss", "pea soup; fp", "madre curry powder hot", "highland barley wine",
    "salad", "rotten eggs", "steamed conch", "rice wine rice syrup", "garam masala",
    "water kefir", "palm sap", "peanut butter", "ground turkey", "sick cider", "onion",
    "austrian hard cheese rind", "litopenaeus vannamei purchased at the supermarket",
    "chinese pickle", "cheese", "raw mutton", "nuruk; korean traditional beverage starter",
    "crab marinated in soy sauce", "pork product", "chinese sauerkraut", "yoghourt",
    "honey", "ham", "bobby veal steak", "10' wieners", "natto", "yogurt", "arugula",
    "gochujang", "cheonggukjang", "cheese starter culture", "corn steep liquor",
    "traditional greek kasseri cheese", "pork & cabbage dumplings", "lettuce", "fish balls",
    "new zealand cheese", "wine", "pork casserole; fp", "the red brine of salted laminaria",
    "beer", "tibetan kefir", "carrot juice", "pickle", "butter starter", "home-made vinegar",
    "product-eggs-raw-whole", "tomato pulp", "frozen whole tilapia", "fish ball", "tomatoes",
    "raw sausage", "beer contaminant", "scoby from kombucha tea", "calf liver", 
    "mustard pickles", "organic spinach", "smoked fish", "cooked ox blood", "baijiu",
    "doubanjiang", "retail turkey", "fresh produce (lettuce; lollo bionda)", "blue berry",
    "apple cider", "fried eel (fish)", "radish and carrot pickled with rice bran and salt",
    "mac and cheese", "shelled pistachios", "marinated fish product", "vinegar pei",
    "10 weeks old 45+ samso cheese", "salami", "chinese traditional sourdough",
    "mixed salads", "tibet kefir", "sugar thick juice", "traditional dairy products",
    "ice cream", "water containing garlic and cabbage", "doenjang(soybean paste)",
    "green chili", "ground beef", "soy sauce mash", "yoghurt", "rye-bran sourdough",
    "koumiss (fermented mare's milk)", "pilsner beer", "blown cheese", "nuruk",
    "pickled green chili peppers", "fish part (slab) sold as silver carp; muscle",
    "wheat sourdough", "cherry", "pork barbeque", "vinegar", "dried tofu", "raw gallnut honey",
    "sauerkraut", "diced lamb", "salt", "commercial salt", "patagonian red wine",
    "korean traditional alcoholic beverage", "tibetan traditional dairy products",
    "fish meal", "bagged lettuce", "wheat beer", "patagonian merlot wine", "imitation crab",
    "apple juice from cider press", "sourdough", "frozen raw shrimp", "palm brown sugar",
    "isolate obtained from campbell's Soup", "brine of stinky tofu", "mead", "sikhae",
    "french dry-type pork sausage", "sichuan pickle vegetables", "emmental",
    "cheese production", "pork dumplings", "home-made koumiss", "korean kefir",
    "dessicated chinese egg powder", "stinky tofu", "leaf vegetable", "retailed chichen",
    "1 month-old fish sauce mash", "10 weeks old samso 45+ cheese", "wiener sausage",
    "doenjang", "flour", "unsweetened puffed rice cereal", "pinot noir wine",
    "alfalfa sprouts", "ground turkey sausage", "turkey; ground", "sponge cake"
    "jogaejeotgal; a traditional korean", "korean soybean paste", "cabbage cv. wk882",
    "retail chichen", "ground red chili pepper", "makgeollli", "cucumbers", "yellowfin tuna scrape",
    "raw milk cheese", "maotai daqu", "jeotgal", "zha-chili", "heat treated raw milk",
    "dairy product", "isolate obtained from campbell's soup", "retail pork", "egg raw whole",
    "frozen yellowfin tuna nakaochi scrape", "product egg raw white", "homemade yoghurt",
    "traditional homemade dairy product", "gochujang(korean red pepper paste)",
    "frozen yellowfin tuna steak", "retail chicken", "retail chichen", "dried roasted pistachio",
    "chicken carcass", "malted barley", "bee honey", "anaheim pepper", "thick broad-bean sauce",
    "chicken product", "pork steak", "algal food product", "ganjang(korean soy sauce)",
    "bakery environment concentrated whipped topping", "congee", "whipping cream",
    "environmental swab sponge bakery mixing bowl outside", "barbecue pork",
    "bakery environment mixing bowl outside", "ground sesame seed paste tahini",
    "fish part (head) sold as asian bighead carp; muscle", "traditional yoghurt",
    "plant derived food stuff; onion; allium cepa", "the brine of stinky tofu",
    "cold-stored modified atmosphere packaged broiler filet strips with the first signs of spoilage"]


livestock_hosts = [
    "companion dog", "fattening pig", "healthy pig", "gallus gallus (chicken)", "ovine (sheep)",
    "canis lupus familiaris breed german shepherd (lovey)", "horse 2/1 (3)", "felis silvestris catus",
    "bos taurus linnaeus", "bos taurus (cow)", "ovine", "shepp", "calves", "dairy cattle", "guinea pig",
    "sus scrofa", "suxs scrofa domesticus", "bos taurus", "piglet", "pig", "cattle", "swiss mouse",
    "cows", "broiler chicken", "cow rumen", "bos taurus (bovine)", "horse", "horse 2/1 (4)",
    "sheep placental tissue", "mouse-c57bl/6j", "swine", "chicken", "goose", "duck", "mammal (dog)",
    "bovine", "capra hircus", "pet dog", "felis catus", "healthy carrier pig", "donkey",
    "gallus gallus", "gallus gallus domesticus", "porcine", "cow", "cat", "equus caballus",
    "mare", "calf", "gallus", "canicola", "hen", "meat duck", "diseased pig", "buffalo",
    "equus caballus ferus", "pork", "sus", "bubalus bubalis", "water buffalo", "dairy herd",
    "goat", "sheep", "equus ferus caballus", "bos bovis", "turkey", "rabbit", "alpaca",
    "bos primigenius taurus", "canine", "laoshan milk goat", "black pig", "domestic animals",
    "canis lupus familiaris", "bos taurus coreanae", "yak", "dog", "equine", "ovis aries",
    "young chicken", "gallus gallus domesticus isa15", "beef cattle", "korean short hair cat",
    "dairy cow", "poultry", "piglets", "lama glama", "ovis aries (domestic sheep)",
    "duck with tremor", "canis lupus familiaris breed parson russell terrier and chihuahua mix",
    "sus domesticus", "canis familiaris", "cattle and buffalo", "bovine (dairy herd)",
    "bird (chicken)", "bird (turkey)", "bos indicus x", "ruminants (bovine)",
    "capra aegagrus hircus (domesticated goat)"]


livestock_isolation_sources = [
    "swine preevisceration carcass swab", "air sacs", "sus scrofa domesticus", "rumen contents",
    "liquid joint sample of pig", "air sac of a pipped embryo", "weaned piglets", "guinea pig",
    "apoblema of swine", "isolated from cattle with blackleg", "psittacosis outbreak", "ropy milk",
    "cattle carcass", "nasal swab taken from a healthy thoroughbred racehorse", "turkey",
    "liquid joint sample of a pig", "sheep", "sheep isolate", "cow; fecal", "faeces of piglets",
    "canine skin", "canine feces", "caecum", "rumen of sheep", "animal - porcine-feces",
    "pig; intestine", "duckling with tremor", "manure", "horse manure", "rumen",
    "choanal cleft of a commercial layer hen with respiratory disease", "pig; small intestine",
    "cattle", "chicken", "farm", "turkey feces", "wastewater from pig manure", "pig leung",
    "chicken feces", "nasal cavities of a calf", "bos taurus", "nasal cavities of a pig",
    "pig", "swine nasal swab", "the intestine membrane of a diarrheic piglet",
    "intestine membrane of a diarrheic piglet", "pig manure", "canine oral cavity",
    "isolated from the intestines of sick birds in the farm", "poultry establishment",
    "slaughtered pig", "poultry", "cat", "duck", "broiler", "pig faeces", "poultry environment",
    "fecal samples from slaughtered sheep", "broiler chicken", "chicken cecal content",
    "dog with mastitis", "goose anus swab", "parent strain cv601 collected from dairy manure",
    "porcine", "cattle hide", "goat; fecal", "rumen fluid", "pig environment", "porcine feces",
    "cow raw milk", "cloacal swab", "skin of a pig", "cattle faeces", "porcine rectal swab",
    "lesion site (lung) of a dead turkey with colibacillosis", "laying hen withcolibacillosis",
    "goose faeces", "cattle slurry", "carcass", "pig fecal", "lamb; fecal", "disease duck",
    "the respiratory tract of a pig with swine respiratory disease", "cattle feces",
    "wool from pakistan", "chicken skin", "raw milk", "swine excrement", "fecal sample of cattle",
    "a nose swab sample of swine origin", "healthy broiler chicken", "chicken manure",
    "alpaca; fecal", "horse with strangles", "fresh raw milk", "sus scrofa: cecal content",
    "swine's gut", "swine final chilled carcass", "cow rumen", "cows", "dog",
    "bovine rumen", "cow dung", "poultry carcass", "air of cow shed", "pig feces",
    "pig feed from feed plant", "fjerkrae", "liver of poult", "horse feces", "swine feces",
    "tissue and/or biological fluid swine", "swine cecum", "raw cow milk", "turkey with septicemic infection",
    "heifer vaginal mucus", "cloaca", "buffalo calf", "aborted piglet fetus",
    "calf", "animal feed", "rectal swab of cattle in slaughterhouse", "horse faces",
    "retail chicken", "canine oral plaque", "bovine pre-evisceration carcass at",
    "novine (bobby calf; fecal)", "a feces sample of chicken origin", "gallus gallus",
    "lung of aborted horse fetus", "beef liver", "swine", "chicken intestine",
    "broiler chick cecum", "isolated from the intestines of sick", "chicken cecum", "bulk pig ear",
    "bulk pig ears", "cattle slurry", "cow feces", "animal-swine-roaster swine",
    "ncsu equine educational unit", "cowshed of a farm", "a feces sample of swine origin",
    "sheep fecal sample", "digested slurry of dairy manure", "diseased pig",
    "a anal swab sample of swine origin", "ear swab from dog", "goat", "healthy weaning piglets",
    "lesion site (lung) of a dead turkey", "chicken dung", "chicken trachea",
    "pooled sheep faecal samples collected from floor of farm", "poultry litter",
    "pooled pig faecal samples collected from floor of farm", "bull calf intestinal microflora",
    "pooled cattle faecal samples collected from floor of farm", "ovine (sheep)",
    "an australian field isolate recovered from choanal cleft of a commercial layer hen with respiratory disease",
    "abomasum content of an aborted sheep's fetus", "cow manure",
    "pig with exfoliative dermatitis",
    "rectal fecal grab samples from a commercial feedlot", "chicken tissues"]


anthropogenic_hosts = ["sewage", "soil (polluted)", "dsmz strain", "laboratory", "lab",
                       "cell culture collection", "mine", "bacillus thuringiensis subsp",
                       "zymobiomics microbial community standard strain"]


anthropogenic_isolation_sources = [
    "artp mutation", "the leaf of flue-cured tobacco", "interstitial fluid", "in vitro substrate",
    "parking lot island in grocery store parking lot", "disinfectant bottle", "kaist",
    "water from shower head", "cell culture", "primary hospital lab", "atcc 35296",
    "hospital universitari germans trias i pujol (badalona)", "digestate", "enriched consortium",
    "hospital universitario son de espases (palma de mallorca)", "elanco animal health",
    "shoe sole from worker", "phosphate mine", "air conditioner condensate drain pipe",
    "re-isolation from atcc 39073", "airborne", "spontaneous mutant", "lab mutant",
    "starvation in luria bertani broth", "thermophilic stage of composting process",
    "dairy farm; merseyside; uk", "flange fragments", "lab mutation of atcc 17978",
    "derived from cip 106327 (collection de l'institute pasteur; paris; france)",
    "strain ev g. gerard; g robic (vaccine strain ev line niieg is mutant strain; was isolated from other stain ev gerard and robic)",
    "thermophilic microbial from zc4 compost from a compost operation in the sao paulo zoo (brazil)",
    "environmental swab", "chimney", "medicine 'bioflor'", "bei isolate", "biomedical source",
    "spontaneous rifampin-resistant isolate from atcc 43816", "gas production well settling pond",
    "cryoprecipitate; transfusion reaction", "atcc reference strain", "alkalinize citrate",
    "disintigrated concrete from an outflow sewer", "cathodes of microbial fuel cell",
    "commercial digestive syrup", "olive production company", "sewage", "chick paper",
    "plastic surface in contract manufacturing organization", "residental environment",
    "slag heap at uranium mining site", "sewage", "pyrite-leaching pilot plant",
    "public walkway near shopping center; kumomi; japan", "equipment", "granule",
    "residential yard with dogs and pecan tree", "distribution center", "biofilm",
    "mining area near to la esperanza in murcia; spain", "vaccine", "biogas plant",
    "superficial sediment of polluted river", "sufu sample", "polyethylene", "filter",
    "dust collector", "shrimp culture pond", "ts-11 based vaccine purchased from merial select",
    "sugar beet juice from extraction installations", "crap-stick", "bioswale",
    "shrimp culture ponds", "lab selection", "dairy culture", "culture", "laboratory",
    "lab-derived ems mutagenesis", "trout farm", "rusted steel wire rope",
    "environmental swab sponge bakery light switches on nw wall", "winery environment",
    "anthropogenic environmental material", "environmental swab sponge bakery hallway",
    "chinese petroleum corporation 77 kaoshiung refinery located in kaoshiung 78 city; taiwan",
    "rusty iron plate", "commercial product serenade (bayer)", "evaporator",
    "enrichment community growing on the high-molecular-weight fraction of a black liquor sample from federal paper board company inc.; augusta; ga", "sunki", "commercial mycovac-l vaccine",
    "sea water aquarium outflow", "vaccinal strain", "retail fish market", "fishmeal",
    "bakery environment - bottom metal shelf on table used for making birthday cakes",
    "gold-copper mining", "commercial strain obtained from novagen", "freshwater aquarium biofilter",
    "petroleum polluted river", "corrosion products from a corroded pipe in an fpso",
    "vinegar factory", "freshwater ras system; circulating water", "traditional dairy",
    "chemical and physical mutagen treatment of atcc 25486", "lactic starter", "starter culture",
    "snottite biofilm from former pyrite mine", "dairy fan", "dsmz", "dsmz isolate",
    "cleaning system aquaculture", "clam larvae aquaculture", "mineral salt medium",
    "combined sewer", "oil", "roundup", "culture mutant", "turbot fish farm",
    "commercial lactic starter", "beet pulp; sugar refinery", "mural paintings",
    "high concentration of fluoride", "fish larvae aquaculture", "sewage water",
    "combined sewer effluent", "feed additive imported from china", "wolfram mine tailing",
    "swab from a follow-up assessment of the tap handles and sink edges after a first disinfection attempt",
    "compost windrow", "residential environment", "esa spacecraft assembly clean room",
    "storage tank", "zc4 compost from a compost operation of sao paulo zoo",
    "swab from a hand-washing sink as part of the hospital routine surveillance program",
    "potash salt dump", "bakery environment - hallway", "r2a medium", "polystyrene",
    "in the nicotine environment", "outdoor built environment", "farm kitchen",
    "clone isolated from the evolution experiment described in ketola et al. 2013 (doi: 10.1111/evo.12148)",
    "the sample was isolated from cell culture after 9 month of microevolution experiment",
    "salt mine", "bioreactor", "sewage & soil", "tung meal", "cooling tower",
    "surface of a polyethylene microplastic particle present in tank 6 of a marine aquarium containing stony-coral fragments and water maintained at 26 degree c", "mutagenisis derived",
    "aquarium water", "biofilm sample", "medium", "dairy slurry uk", "landfill",
    "detritus agregates formed in tank 6 of a marine aquarium containing stony-coral fragments and water maintained at 26 degree c", "duedonoscope instrument", "air in a school dining room",
    "non-filtered water from the water column of tank 6 of a marine aquarium containing stony-coral fragments. water maintained at 26 degree c", "cheddar cheese factory", "mine", "brewing yeast sample",
    "anabaena culture", "dairy environment", "genome shuffling", "jugular catheter",
    "telephone of nurse station", "shower 3", "geothermal power plant", "sink aerator",
    "shelf", "brewery environment", "lab strain", "coffee cup", "culture collection",
    "air in laboratory", "sp4 medium; universite de rennes", "industrial", "knife", "wet market",
    "hospital", "open pond on an algae farm", "carbonate chimney", "concrete",
    "sink handle in icu room in military hospital", "hospital environment",
    "sf9 cell culture media", "laboratory", "mold-colonized wall of an indoor",
    "power plant biotrickling filter", "biopesticide dipel df", "biopesticide delfin",
    "biopesticide novodor 3fc", "biopesticide solbac", "biopesticide xentari",
    "alcohol foam dispenser in hospital intensive care unit", "oilfield", "book surface",
    "bedside light switch in hospital intensive care unit", "1970 vaccine stock",
    "shower 2", "room 7", "research laboratory", "mycoplasma culture contaminant",
    "anaerobic digestion reactor", "anthropogenic", "cattle slaughter plant",
    "heated-cooler unit water tank", "beer keg", "woodchip bioreactor", "environmental sponge",
    "wastewater treatment plant effluent", "ac condensate", "industry",
    "bed sheets", "automobile air-conditioning evaporator", "agricultural settling lagoon",
    "biocathode mcl", "envo:01000905", "indoor air", "cultured in erlenmeyer flask",
    "underground farm slurry reservoir", "thermal field", "crude oil", "lb (luria broth; bd; usa)",
    "sanger centre via imperial college", "market", "maiket", "cutting board",
    "pulp of a gold-containing sulfide concentrate", "zebrafish tank detritus",
    "polycyclic aromatic hydrocarbon", "agar plate", "brewery", "latex", "brewery-associated surface",
    "bakery environment - assembly production room", "hot-gas well; coatings inside tube",
    "anodic biofilm of glucose-fed microbial", "milk powder production facility",
    "coke and gas plant treatment facilities", "liquid-air interface biofilm",
    "rubber production plant territory", "laboratory stock", "air conditioning system",
    "ensuite 7/8", "river sediment polluted by acid mine drainage",
    "environmental swab sponge bakery assembly production room",
    "microbial mats from zloty stok gold mine", "acid mine decant and tailings from uranium mine",
    "sedimentation pond in a zinc factory", "acid mine drainage", "filter from dairy farm",
    "orfrc groundwater during biostimulation for uranium bioreduction",
    "sink drain", "chemical manufacturing sites", "borehole hdn1; spa", "biofilm boat",
    "e-waste recycling site", "paper mill kaolin", "municipal wwtp",
    "np30 nitrate/perchlorate-reducing mixed", "dust collector of pigpen",
    "acid mining effluent decantation pond", "fish farm", "subarctic landfill",
    "anaerobic digester", "cooling tower water", "dairy", "uninoculated hep-2 tissue culture",
    "tattoo ink", "environment of small animal veterinary clinic",
    "derivative of strain mm294 except the mutated alleles reca1 and gyra96",
    "isolated from commercially available", "contaminated wipes", "protoplast breeding",
    "arjo bathroom", "paper surface", "composted garbage", "woad vat",
    "floor surface of biological laboratory", "experimental exposure to phage",
    "hot water tap; geest office building", "industrial building air scrubber liquid",
    "iron water pipe", "sterile tools", "dental water pipeline", "sungai pinang; penang; malaysia",
    "tar at the shipwreck site of tanker haven", "human septic tank", "monitor panel",
    "western north pacific station s1", "ventilator", "coal-cleaning residues", "cured hides",
    "fresh-cut produce processing plant", "burns unit surveillance", "coal mine",
    "chilean kraft-pulp mill effluents", "automobile air conditional evaporator",
    "hospital sink", "oil-immersed sample from guaymas basin", "oil contaninated soil",
    "underground horizons of a flooded mine in russia", "canadian salted buffalo hide",
    "red heat in salted hides (spoiled fish)", "colorado serum co.; anthrax spore vaccine",
    "electroporation of y pestis kim6+ in pcd1ap+", "solar panel array [envo:01000867]",
    "denitrifying; sulfide-oxidizing effluent-treatment plant", "vaccine isolate",
    "high motility on swim plates; streptomycin resistant", "insecticide factory",
    "enrichment cultures from ucc lynda-stan", "from type material", "drain", "takara",
    "base of single tree in sidewalk along rollins rd. evidence of many dogs; boston; mass",
    "monkey kidney tissue-culture fluids of the fh strain (eaton agent virus)",
    "conjugation assay", "lab", "culture maintained in leon; missing plasmid pscl3",
    "probitic products", "washroom sink in hospital intensive care unit",
    "doe field research center at oakridge; tennessee; area 5; well fw507",
    "hygromycin b antibiotic bottle", "paper pulp mill", "residential yard",
    "biofilm reactor", "nutrient broth"
]
sediment_hosts = ["ocean sediments", "sediment"]

sediment_isolation_sources = [
    "mud", "black mud", "hydrothermal vent area derived sediment", "envo:01001050", "tidal flat",
    "marine mud", "muddy water", "pit mud", "tidal marsh", "freshwater mud", "a tidal flat",
    "mangrove swamp", "sea mud", "mud flat", "the pit mud of a chinese liquor", "sea-tidal flat",
    "tidal flat sample", "pink berry consortia", "submerged sand bank", "ocean sediments",
    "muddy water", "organically rich wetland", "bog", "river clay", "mangrove habitat samples",
    "mangrove wetland ecosystem", "subseafloor basaltic crust", "lake mainaki salt",
    "mud from a salt lagoon", "envo:00002113", "envo:00002007; envo:00002007"]

human_hosts = ["ethnic koreans living in china", "young woman","nist mixed microbial rm strain",
               "infant", "soldier"]

human_isolation_sources = [
    "isolate from baby feces", "bile", "faecal swabs", "healthy gingiva", "cecal contents", "ucc isolate",
    "ventricular fluid", "body surface", "skin sampling", "bile bronchoalveolar lavage fluid", "carious dentine",
    "cephalorachidian fluid", "adult feces", "gingival margin", "oral swab", "respiratory sample", "knee fluid",
    "microbial feature", "abdominal skin swab", "amniotic fluid", "subcutaneous granuloma",
    "intestine of adult", "intestinal tract", "intestinal contents", "oviduct", "cutaneous", "vaginal swab",
    "granulomatous lesion", "fecal sample of a healthy adult", "aerosol sample", "ear", "nares",
    "bronchoscopy washings", "diarrhoea", "tracheal suction", "corneal ulcer", "oral cavity",
    "cerebrospinal fluid", "broncho-alveolar lavage", "cystic fibrosis lung", "scabs", "faecal sample",
    "pediatric cohort consisted of 67 unique buccal molar plaque samples", "leg ulcer", "gingivitis",
    "genitourinary tract", "new guinea indigenes", "leg", "cf sputum", "spleen", "ear infection",
    "gastric antrum from amerindian resident", "faeces / diarrhoea", "centenarian fecal sample",
    "massachusetts listeriosis outbreak", "pharyngeal mucosa", "pre-capsular lymph node",
    "bronchoalveolar lavage", "cerivx", "vagina", "tonsil", "carious teeth", "nose", "pharyngeal tonsil",
    "baby feces", "high vaginal swab", "vagina of healthy woman", "thumb", "abscess", "bone",
    "from patient with pneumonia", "blood", "stool", "saliva", "skin swab", "homo sapiens",
    "wound", "dental abscess", "compound fracture", "dog bite", "oral mucosa", "trachea",
    "rectal swab", "skin", "fecal", "oral", "nosocomial environment", "gut", "pus", "biopsy",
    "human tissue biopsy", "korean adult feces", "groin", "infant's throat", "post nasal swab",
    "sputum (cystic fibrosis patient)", "diarrheal patient", "chest pus", "deep liver", "lung tissue",
    "clinical patient", "intestine", "respiratory tract", "human clinical sample", "lymph node",
    "vomit of food poisoning patient", "faecal", "respiratory", "materia alba", "nasopharynx",
    "human skin swab", "vomit from a food poisoning case", "superficial wound", "conjuctival",
    "china: the first affiliated hospital of zhejiang university", "sputum", "leg wound",
    "human fecal sample", "human blood", "blood sample", "urinary tract infection", "urogenital",
    "bronchial lavage", "human feces", "anus swab", "feces; fp", "nasopharyngeal mucus",
    "wound infection", "human", "stool of a patient presenting with", "breast abscess",
    "gastric", "infant feces", "eye", "semen", "inguinal/rectal swab", "foot wound",
    "human feces (woman; 24 years old)", "human stool", "urine", "faeces", "inflamed gingiva",
    "skin abscess", "hepatic abscess", "burn wound", "medical specimen", "subgingival plaque",
    "tibia/osteomyelitis", "patient", "lung sample", "infant faecal", "finger", "dental caries",
    "isolated from blood of a patient with", "16-year-old boy with bubonic plague",
    "fecal sample", "tissue", "ascites", "vaginal secretions", "diarrhea", "genital tract",
    "kidney", "liver", "nose swab", "throat swab", "mouth", "dental plaque", "gastric antrum",
    "blood of a hospitalized patient", "brain", "human feces (woman; 60 years old)",
    "foodborne disease surveillance", "supragingival dental plaque", "gingival sulcus",
    "korean infant feces", "healthy adult male gut", "oropharynx", "human foot", "vein blood",
    "perirectal swab", "rectal", "anal swab", "human skin tissue", "a diarrheal child",
    "blood specimen", "rectal fecal grab samples from a", "throat", "cystic fibrosis patient",
    "hospitalized patients", "from patient with wound infection", "human gi tract",
    "saliva of a healthy person",
    "human case of meningitis", "supragingival plaque; periodontitis", "spinal fluid",
    "cancer patient stool", "milker's hand", "hip replacement", "oropharnynges", "infected wound",
    "biopsy of antral stomach region from european patient with peptic ucler and chronic gastritis disease"]

na_hosts = ["environmental", "jamiecosley", "environment", "not available: not collected", "microbial",
            "ucc strain", "natural / free-living", "obscured"]

na_isolation_sources = [
    "environmental material [envo:00010483]", "japan", "mssing", "nanchang", "environmetal",
    "envo:00010483", "not provided; submitted under migs 2.1", "huangshui", "shanghai", "ucc strain",
    "enviornmental", "host's whole body", "bacterial consortium", "missing", "unknown", "environmental",
    "afb-diseased colony", "environment swab", "physical", "free-living", "cell culture", "bacteria",
    "biological fluid", "not available: to be reported later", "environment", "unavailable",
    "r.j.roberts", "whole body", "atcc isolate", "feed", "not available: not collected", "tianjin",
    "isolate obtained from atcc:7955 nca", "kimoto", "ncimb", "valley", "conjugation", "kyoto",
    "obscured", "peter dedon mit", "swab", "dsm:3754", "laboratory isolate", "tsoundzou", "surface",
    "environmental surface", "ghana", "minnesota", "jiangsu", "not available", "iam12617", "singapore", 
    "usa: wa: seattle", "not available: to be reported later", "not known", "bone powder",
    "greensboro; alabama", "n.a.", "cameron currie (currie@bact.wisc.edu)", "no",
    "terry c. hazen (tchazen@lbl.gov)",
    "derek lovley (dlovley@microbio.umass.edu)", "atcc strain", "atcc", "china: tianjin",
    "laboratory strain derived from cip", "local geographical source", "neb433", "not isolated",
    "swab with brown-gray powder", "epidemic", "anaerobic environments", "environmental", "neb269",
    "environmental swabs", "uruguay", "coconi", "csf", "balf", "china: zibo; shandong",
    "integrative microbiology research center; south china agricultural university (scau); tianhe district; guangzhou china",
    "gram-negative and facultative anaerobes; with flagellated and non-spore-forming",
    "wild type strain isolated from a natural source", "beijing", "poland; warsaw area",
    "bal", "unkown", "other", "to be completed", "agr", "enviornmental"]

#######################################################################################
#Script Begins Here
#######################################################################################


print("Annotation_Accession,host,isolation_source,Annotation")  ## print header

with open("../results/gbk-annotation-table.csv", "r") as annotation_fh:
    for i, line in enumerate(annotation_fh):
        if i == 0: continue ## skip the header.
        line = line.strip()  ## remove trailing newline characters
        fields = line.split(',')  ## split line into an array
        annotation_accession, original_host, original_isolation_source = fields
        ## turn host and isolation source to lower case (NA goes to na!)
        host = original_host.lower()
        isolation_source = original_isolation_source.lower()
        annotation = "blank" ## default value, so we know whether it is set later on.
        
        ## NA
        if host == "na" and isolation_source == "na":
            annotation = "NA" ## the output NA should be upper-case.
        if host == "na" and isolation_source in na_isolation_sources:
            annotation = "NA" ## the output NA should be upper-case.
        elif isolation_source in na_isolation_sources and host == "na":
            annotation = "NA"
        elif host in na_hosts and isolation_source == "na":
            annotation = "NA"
        elif host in na_hosts and isolation_source in na_isolation_sources:
            annotation = "NA"

        ## Terrestrial
        elif host == "na" and isolation_source in terrestrial_isolation_sources:
            annotation = "Terrestrial"
        elif "volcano" in isolation_source:
            annotation = "Terrestrial"
        elif "acidic" in isolation_source and "wetland" not in isolation_source:
            annotation = "Terrestrial"
        elif "solfataric" in isolation_source:
            annotation = "Terrestrial"
        elif "subsurface" in isolation_source and "sediment" not in isolation_source:
            annotation = "Terrestrial"
        elif "rifle well" in isolation_source:
            annotation = "Terrestrial"
        elif "soda lake" in isolation_source:
            annotation = "Terrestrial"
        elif "hypersaline" in isolation_source:
            annotation = "Terrestrial"
        elif "hot spring" in isolation_source and "soil" not in isolation_source:
            annotation = "Terrestrial"
        elif host == "air":
            annotation = "Terrestrial"

        ## Soil
        elif "soil" in isolation_source and "contamina" not in isolation_source and "pollute" not in isolation_source and "contaninated" not in isolation_source and "sewage" not in isolation_source:
            annotation = "Soil"
        elif host == "na" and isolation_source in soil_isolation_sources:
            annotation = "Soil"
        elif "soil" in host and "contamina" not in host and "pollute" not in host:
            annotation = "Soil"
        elif "rhizosphere" in isolation_source or "rhizosphere" in host:
            annotation = "Soil"
        elif host in soil_hosts:
            annotation = "Soil"

        ## Agriculture
        elif host in agri_hosts and host not in plant_hosts and isolation_source not in food_isolation_sources:
            annotation = "Agriculture"
        elif isolation_source in agri_isolation_sources and host not in animal_hosts:
            annotation = "Agriculture"
        elif "sativa" in host:
            annotation = "Agriculture"
        elif "citrus" in host:
            annotation = "Agriculture"
        elif "solanum" in host:
            annotation = "Agriculture"
        elif "tomato" in host:
            annotation = "Agriculture"

        ## Plant-host
        elif host in plant_hosts:
            annotation = "Plants"
        elif isolation_source in plant_isolation_sources:
            annotation = "Plants"

            
        ## Marine
        elif host in marine_hosts:
            annotation = "Marine"
        elif isolation_source in marine_isolation_sources:
            annotation = "Marine"
        elif host == "na" and " sea " in isolation_source and "sediment" not in isolation_source and "mud" not in isolation_source and "yogurt" not in isolation_source:
            annotation = "Marine"
        elif "seawater" in isolation_source and host == "na" and "aquaculture" not in isolation_source:
            annotation = "Marine"
        elif "hydrothermal" in isolation_source and host == "na" and "sediment" not in isolation_source:
            annotation = "Marine"
        elif ("marine" in isolation_source and "sediment" not in isolation_source and "mud" not in isolation_source
              and host not in animal_hosts and isolation_source not in terrestrial_isolation_sources
              and "sponge" not in isolation_source and "aquarium" not in isolation_source):
            annotation = "Marine"

            
        ## Freshwater
        elif host in freshwater_hosts:
            annotation = "Freshwater"
        elif host == "rhodobacter sphaeroides 2.4.1":
            annotation = "Freshwater"
        elif host == "cereibacter sphaeroides 2.4.1":
            annotation = "Freshwater"
        elif host == "environment" and isolation_source == "water":
            annotation = "Freshwater"
        elif host == "environmental" and isolation_source == "water":
            annotation = "Freshwater"
        elif ("water" in isolation_source and isolation_source not in anthropogenic_isolation_sources
              and isolation_source not in sediment_isolation_sources
              and "sea" not in isolation_source and "waste" not in isolation_source and "oil" not in isolation_source and "aquarium" not in isolation_source
              and host == "na" and "ocean" not in isolation_source and "sediment" not in isolation_source
              and isolation_source not in food_isolation_sources and "sal" not in isolation_source and "kimchi" not in isolation_source):
            annotation = "Freshwater"
        elif "spring" in isolation_source and "oil" not in isolation_source and "sediment" not in isolation_source:
            annotation = "Freshwater"
        elif ("lake" in isolation_source and "sediment" not in isolation_source and "sal" not in isolation_source
              and "soda" not in isolation_source and "oil" not in isolation_source and "silt" not in isolation_source):
            annotation = "Freshwater"
        elif isolation_source in freshwater_isolation_sources:
            annotation = "Freshwater"

            
        ## Food
        elif host in food_hosts:
            annotation = "Food"
        elif "pickle" in host:
            annotation = "Food"
        elif isolation_source in food_isolation_sources:
            annotation = "Food"
        elif (host == "na" and "milk" in isolation_source and "hand" not in isolation_source
              and "raw" not in isolation_source and "facility" not in isolation_source):
            annotation = "Food"
        elif host == "na" and "food" in isolation_source and "vomit" not in isolation_source and "sludge" not in isolation_source:
            annotation = "Food"
        elif host == "na" and "beef" in isolation_source:
            annotation = "Food"
        elif "cheese" in isolation_source:
            annotation = "Food"
        elif "pickle" in isolation_source:
            annotation = "Food"
        elif "kimchi" in isolation_source:
            annotation = "Food"
        elif "meat" in isolation_source and "facility" not in isolation_source:
            annotation = "Food"
        elif "ferment" in isolation_source:
            annotation = "Food"
        elif host == "na" and "chicken" in isolation_source and isolation_source not in livestock_isolation_sources:
            annotation = "Food"
        elif host == "na" and "yogurt" in isolation_source:
            annotation = "Food"

            
        ## Livestock
        elif host in livestock_hosts and isolation_source not in food_isolation_sources:
            annotation = "Livestock"
        elif isolation_source in livestock_isolation_sources and host not in animal_hosts:
            annotation = "Livestock"
        elif "bovine" in isolation_source:
            annotation = "Livestock"
        elif "scrofa" in host:
            annotation = "Livestock"
        elif "placental tissue" in isolation_source:
            annotation = "Livestock"
        elif isolation_source == "turkey" and host == "na":
            annotation = "Livestock"

            
        ## Animal-host
        elif host in animal_hosts:
            annotation = "Animals"
        elif isolation_source in animal_isolation_sources:
            annotation = "Animals"
        elif "bemisia" in host: ## as in bemisia tabaci.
            annotation = "Animals"


        ## Anthropogenic-environment
        elif host in anthropogenic_hosts:
            annotation = "Human-impacted" 
        elif host == "na" and isolation_source in anthropogenic_isolation_sources:
            annotation = "Human-impacted"
        elif host == "environment" and isolation_source in anthropogenic_isolation_sources:
            annotation = "Human-impacted"
        elif host == "environmental" and isolation_source in anthropogenic_isolation_sources:
            annotation = "Human-impacted"
        elif host == "sediment" and isolation_source in anthropogenic_isolation_sources:
            annotation = "Human-impacted"
        elif host == "na" and "aquaculture" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "waste" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "sewage" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "sludge" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "bakery environment" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "modified" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "contamina" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "digester" in isolation_source:
            annotation = "Human-impacted"
        elif "oil" in isolation_source:
            annotation = "Human-impacted"
        elif "cell line" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "industr" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "construct" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "artificial" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "facility" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "engineer" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "culture" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "agar" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "reactor" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "laboratory" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "experiment" in isolation_source:
            annotation = "Human-impacted"
        elif host == "na" and "biopesticide" in isolation_source:
            annotation = "Human-impacted"
        elif "polluted soil" in isolation_source:
            annotation = "Human-impacted"
        elif host == "gram-positive bacteria" or host == "bacteria" or host == "sludge":
            annotation = "Human-impacted"
        elif "probiotic" in isolation_source:
            annotation = "Human-impacted"
        elif host ==  "sofa":
            annotation = "Human-impacted"

            
        ## Sediment
        elif "sediment" in isolation_source and "contamina" not in isolation_source and "pollute" not in isolation_source:
            annotation = "Sediment"
        elif host in sediment_hosts:
            annotation = "Sediment"
        elif host == "na" and isolation_source in sediment_isolation_sources:
            annotation = "Sediment"
        elif host == "marine sediment":
            annotation = "Sediment"
        elif "mud" in isolation_source:
            annotation = "Sediment"
        elif "tidal flat" in isolation_source:
            annotation = "Sediment"
        elif "silt" in isolation_source:
            annotation = "Sediment"
        elif "wetland" in isolation_source:
            annotation = "Sediment"
        elif "estuary" in isolation_source:
            annotation = "Sediment"
            

        ## Human-host
        elif "sapiens" in host:
            annotation = "Humans"
        elif host != "na" and host not in livestock_hosts and host not in animal_hosts and "homo" in host:
            annotation = "Humans"
        elif host != "na" and host not in livestock_hosts and host not in animal_hosts and "human" in host:
            annotation = "Humans"
        elif host != "na" and host not in livestock_hosts and host not in animal_hosts and host == "human":
            annotation = "Humans"
        elif host not in livestock_hosts and host not in animal_hosts and isolation_source in human_isolation_sources:
            annotation = "Humans"
        elif host in human_hosts:
            annotation = "Humans"
        elif host == "na" and isolation_source == "feces":
            annotation = "Humans"
        elif host == "na" and "blood" in isolation_source:
            annotation = "Humans"
        elif host == "na" and "lymph" in isolation_source:
            annotation = "Humans"
        elif host == "na" and "abscess" in isolation_source:
            annotation = "Humans"
        elif host == "na" and "patient" in isolation_source:
            annotation = "Humans"
        elif host == "na" and "infant" in isolation_source:
            annotation = "Humans"
        elif host == "na" and "human" in isolation_source:
            annotation = "Humans"
        elif host == "na" and "stool" in isolation_source:
            annotation = "Humans"
        elif host == "na" and "clinical" in isolation_source:
            annotation = "Humans"
        elif host == "na" and isolation_source == "lung":
            annotation = "Humans"

            
        print(','.join([annotation_accession, original_host, original_isolation_source, annotation]))

