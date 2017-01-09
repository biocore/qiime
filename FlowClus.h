/*
  John M. Gaspar (jsh58@unh.edu)
  June 2013

  Header for FlowClus.c
*/

// constants
#define MAX_SIZE   65528  // maximum length (characters) for each line of input
#define HEADER     20     // maximum header length
#define NUC        4      // number of nucleotides (ACGT)
#define MIDPRIM    100    // maximum mid tag - primer length
#define DELIM      " \n"  // delimiter for flow values
#define CSV        ",\t"  // delimiter for primer - mid tag file
#define SEP        "-"    // separator for header information
#define PER        "_"    // separator for QIIME, Perseus output
#define COM        ","    // separator for chimera mapping file
#define END        -1.0f  // tag to indicate end of flowgram
#define MIN        0.50f  // maximum flow value not to call a base (rounding down)
                          //   -- change to 0.49f to round up

// label names in sff.txt file
#define NUMFL      "  # of Flows"
#define CHARS      "  Flow Chars"
#define CQL        "  Clip Qual Left"
#define CQR        "  Clip Qual Right"
#define FLOWG      "Flowgram"
#define FLOWI      "Flow Indexes"
#define BASE       "Bases"
#define QUAL       "Quality Scores"
#define COL        ":"

// label names in master file
#define PRIMER     "primer"
#define MIDTAG     "midtag"
#define REVERSE    "reverse"

// default values
#define DEFFLOWEXT ".flow"          // default file extension for cleaned flowgrams
#define DEFFEXT    ".den"           // default file extension for denoised flowgrams
#define DEFMEXT    ".map"           // default file extension for mapping files
#define DEFCHEXT   ".chfasta"       // default file extension for output files for UCHIME
#define DEFCHMAP   ".chmap"         // default file extension for mapping files for UCHIME
#define DEFSDFILE  "stddev.txt"     // default input file containing standard deviations
#define DEFMAXFLOW 19.99f           // default maximum flow value
#define DEFDENFILE "denoised.fasta" // default output file for denoised fasta

/***** command-line options *****/
#define HELP       "-h"  // print usage

// analysis options
#define CLEANOPT   "-a"  // option to clean only (eliminate, truncate)
#define DENOPT     "-b"  // option to denoise only
#define BOTHOPT    "-ab" // option to do both cleaning and denoising (default)

// input/output files
#define MASTERFILE "-m"  // input master file, containing primer and mid tag sequences
#define SFFFILE    "-i"  // input sff.txt file -- required if filtering
#define OUTFILE    "-e"  // output fasta file after filtering
#define DENFASTA   "-o"  // output fasta file after denoising
#define NOMIDOPT   "-x"  // option to produce "QIIME-style" output fasta file(s)
                         //   (no mid tag - primer sequences)
#define FLOWEXT    "-f"  // file extension for filtered flowgrams
#define DENPOPT    "-v"  // option to produce consensus flowgram and mapping files
#define DENFEXT    "-vf" // file extension for denoised flowgrams
#define DENMEXT    "-vm" // file extension for mapping files after denoising
#define ERRFILE    "-c"  // output file for filtering data
#define MISSFILE   "-d"  // output file for misses
#define SDFILE     "-sd" // input file containing standard deviations for each flow value
#define CHIMOPT    "-ch" // option to produce output fasta files for de novo chimera-checking
#define UCHEXT     "-cu" // file extension for UCHIME output fasta files
#define UMAPEXT    "-cm" // file extension for chimera-checking mapping files
#define PEREXT     "-cp" // file extension for Perseus output fasta files

// sequence analysis
#define MIDMIS     "-em" // number of mismatches to mid tag sequence to allow
#define PRIMMIS    "-ep" // number of mismatches to primer sequence to allow
#define MINSLEN    "-l"  // minimum sequence length
#define MAXSLEN    "-L"  // maximum sequence length
#define MAXTRLEN   "-t"  // maximum length at which to truncate a sequence
#define MAXAMBIG   "-N"  // maximum number of ambiguous bases
#define OKAMBIG    "-n"  // number of ambiguous bases to allow before truncating
#define MAXHOMO    "-G"  // maximum homopolymer length to allow
#define OKHOMO     "-g"  // maximum homopolymer length to allow before truncating
#define REVMOPT    "-r"  // option to remove opposite primer if it's found
#define REVQOPT    "-rq" // option to require opposite primer in a read
#define REVMIS     "-er" // number of mismatches to reverse primer sequence to allow

// quality-score analysis
#define AVGQUAL    "-s"  // average quality score
#define WINDOWLEN  "-wl" // length for sliding window of quality scores
#define WINDOWAVG  "-wq" // average quality score for sliding window
#define WINDOWOPT  "-wx" // option to throw out a read if there is a bad window

// flowgram analysis
#define MAXFLOW    "-u"  // maximum flow value -- larger values will be changed to this
#define MINFLEN    "-lf" // minimum flowgram length
#define MAXFLEN    "-Lf" // maximum flowgram length
#define MININT     "-p"  // minimum flow value of interval
#define MAXINT     "-q"  // maximum flow value of interval
#define MAXNVAL    "-z"  // maximum flow value to truncate at
#define NOFLOW     "-y"  // truncate flowgram if this flow value is not
                         //   reached for 4 consecutive flows

// denoising options
#define CINTER     "-j"  // constant value for denoising
#define ZINTER     "-k"  // number of std devs for denoising
#define TRIEOPT    "-tr" // option to denoise using a trie

// error messages
#define ERRMEM     0
#define MERRMEM    "Cannot alloc memory"
#define ERRCLOSE   1
#define MERRCLOSE  "Cannot close file"
#define ERRPARAM   2
#define MERRPARAM  ": cannot find parameter"
#define ERRFLOAT   3
#define MERRFLOAT  ": cannot convert to float"
#define ERRINT     4
#define MERRINT    ": cannot convert to int"
#define ERROPENR   5
#define MERROPENR  ": cannot open file for reading"
#define ERROPENW   6
#define MERROPENW  ": cannot open file for writing"
#define ERRINVAL   7
#define MERRINVAL  ": invalid parameter or usage"
#define ERREXIST   8
#define MERREXIST  ": file already exists"
#define ERRLOAD    9
#define MERRLOAD   "Cannot load indexes or quality scores"
#define ERRNEED    10
#define MERRNEED   ": missing required file"
#define ERRSLEN    11
#define MERRSLEN   "Invalid min/max sequence length"
#define ERRWIND    12
#define MERRWIND   "Invalid quality window values"
#define ERRINTER   13
#define MERRINTER  "Invalid min/max interval flow values"
#define ERRFLEN    14
#define MERRFLEN   "Invalid min/max flowgram length"
#define ERRMAXFL   15
#define MERRMAXFL  "Invalid max noisy flow value"
#define ERRSDVAL   16
#define MERRSDVAL  "Cannot load distance value for denoising"
#define ERRDEN     17
#define MERRDEN    "Must specify either constant value or factor for denoising"
#define ERRDENVAL  18
#define MERRDENVAL "Invalid denoising value"
#define ERRPRIM    19
#define MERRPRIM   "Invalid base in primer"
#define ERRPREP    20
#define MERRPREP   ": cannot repeat primer name"
#define ERRMREP    21
#define MERRMREP   ": cannot repeat midtag name within a primer"
#define ERRNFLOWS  22
#define MERRNFLOWS "Invalid header information in flowgram file"
#define ERRHEAD    23
#define MERRHEAD   "Invalid read header"
#define ERRORDER   24
#define MERRORDER  "Invalid flow order"
#define ERRMID     25
#define MERRMID    ": cannot find mid tag"
#define ERRREV     26
#define MERRREV    "Cannot specify both remove and require reverse primer"
#define ERRNOREV   27
#define MERRNOREV  ": no reverse primer specified"
#define ERRMISM    28
#define MERRMISM   "Invalid number of primer mismatches"
#define ERRMAXF    29
#define MERRMAXF   "Invalid absolute maximum flow value"
#define UNKNOWN    "Unknown error"

// elimination/truncation criteria
#define ETCAT      2     // number of designations -- elim or trunc
#define ELIM       0
#define TRUNC      1
#define SELIM      "Reads eliminated"
#define STRUNC     "Reads truncated"
#define COUNT      "Reads analyzed"
#define MATCH      0
#define PRINT      1
#define SMATCH     "Mid-primer matches"
#define SPRINT     "Reads printed"
#define ERRCAT     16    // number of categories
#define NOERR      0
#define EMINSLEN   1
#define DMINSLEN   "Min. sequence length"
#define EMAXSLEN   2
#define DMAXSLEN   "Max. sequence length for elimination"
#define EMAXTRLEN  3
#define DMAXTRLEN  "Max. sequence length for truncation"
#define EMAXAMBIG  4
#define DMAXAMBIG  "Max. ambiguous bases allowed"
#define EOKAMBIG   5
#define DOKAMBIG   "Max. ambiguous bases allowed before truncation"
#define EMAXHOMO   6
#define DMAXHOMO   "Max. homopolymer length allowed"
#define EOKHOMO    7
#define DOKHOMO    "Max. homopolymer length allowed before truncation"
#define EREVERSE   8
#define DREVERSE   "Reverse primer removed"
#define EAVGQUAL   9
#define DAVGQUAL   "Min. average quality score"
#define EWINDOW    10
#define DWINDOW    "Min. window quality score"
#define EMINFLEN   11
#define DMINFLEN   "Min. flowgram length"
#define EMAXFLEN   12
#define DMAXFLEN   "Max. flowgram length"
#define EFLOWINT   13
#define DFLOWINT   "Noisy flow interval"
#define EMAXNVAL   14
#define DMAXNVAL   "Max. flow value"
#define ENOFLOW    15
#define DNOFLOW    "Four consecutive flows below min."
#define NA         "n/a"
#define TOTAL      "Total\n"

// structs
typedef struct read {
  int length;
  int start;
  float* flow;
  char* header;
  struct midtag* mid;
  struct read* next;
} Read;

typedef struct cluster {
  float* flows;
  int* weight;
  Read* first;
  Read* lon; // longest read
  struct cluster* next;
} Cluster;

typedef struct midtag {
  char* name;
  char* seq;
  int num;
  struct midtag* next;
  struct primer* prim;
} Midtag;

typedef struct node {
  float* flow;
  int st;
  int end;
  int num;
  Read* first;
  struct node* next;
  struct node* child;
} Node;

typedef struct primer {
  char* name;
  char* seq;
  char* rev;
  Midtag* first;
  FILE* out;
  FILE* den;
  FILE* map;
  FILE* per;
  FILE* cmap;
  struct primer* next;
  Cluster* head;
  Cluster* tail;
  Read* dummy;
  Read* prev;
  Node* root;
} Primer;
