from Types cimport *
from libcpp cimport bool
from libcpp.set cimport set as libcpp_set
from String cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/ResidueModification.h>" namespace "OpenMS":
    
    cdef cppclass ResidueModification "OpenMS::ResidueModification":
        # wrap-hash:
        #   getFullId().c_str()

        ResidueModification() nogil except +
        ResidueModification(ResidueModification) nogil except +

        bool operator==(ResidueModification & modification) nogil except +
        bool operator!=(ResidueModification & modification) nogil except +

        void setId(String & id_) nogil except +
        String  getId() nogil except +
        void setFullId(String & full_id) nogil except +
        String  getFullId() nogil except +
        void setUniModAccession(String & id_) nogil except +
        String  getUniModAccession() nogil except +
        void setPSIMODAccession(String & id_) nogil except +
        String  getPSIMODAccession() nogil except +
        void setFullName(String & full_name) nogil except +
        String  getFullName() nogil except +
        void setName(String & name) nogil except +
        String  getName() nogil except +
        void setTermSpecificity(TermSpecificity term_spec) nogil except +
        void setTermSpecificity(String & name) nogil except +
        TermSpecificity getTermSpecificity() nogil except +
        String getTermSpecificityName(TermSpecificity ) nogil except +
        void setOrigin(char origin) nogil except +
        char  getOrigin() nogil except +
        void setSourceClassification(String & classification) nogil except +
        void setSourceClassification(SourceClassification classification) nogil except +
        SourceClassification getSourceClassification() nogil except +
        String getSourceClassificationName(SourceClassification classification) nogil except +
        void setAverageMass(double mass) nogil except +
        double getAverageMass() nogil except +
        void setMonoMass(double mass) nogil except +
        double getMonoMass() nogil except +
        void setDiffAverageMass(double mass) nogil except +
        double getDiffAverageMass() nogil except +
        void setDiffMonoMass(double mass) nogil except +
        double getDiffMonoMass() nogil except +
        void setFormula(String & composition) nogil except +
        String  getFormula() nogil except +
        void setDiffFormula(EmpiricalFormula & diff_formula) nogil except +
        EmpiricalFormula  getDiffFormula() nogil except +
        void setSynonyms(libcpp_set[ String ] & synonyms) nogil except +
        void addSynonym(String & synonym) nogil except +
        libcpp_set[ String ]  getSynonyms() nogil except +
        void setNeutralLossDiffFormula(EmpiricalFormula & loss) nogil except +
        EmpiricalFormula  getNeutralLossDiffFormula() nogil except +
        void setNeutralLossMonoMass(double mono_mass) nogil except +
        double getNeutralLossMonoMass() nogil except +
        void setNeutralLossAverageMass(double average_mass) nogil except +
        double getNeutralLossAverageMass() nogil except +
        bool hasNeutralLoss() nogil except +

cdef extern from "<OpenMS/CHEMISTRY/ResidueModification.h>" namespace "OpenMS::ResidueModification":
    cdef enum TermSpecificity "OpenMS::ResidueModification::TermSpecificity":
        #wrap-attach:
        #    ResidueModification
        ANYWHERE
        C_TERM
        N_TERM
        PROTEIN_C_TERM
        PROTEIN_N_TERM
        NUMBER_OF_TERM_SPECIFICITY

cdef extern from "<OpenMS/CHEMISTRY/ResidueModification.h>" namespace "OpenMS::ResidueModification":
    cdef enum SourceClassification "OpenMS::ResidueModification::SourceClassification":
        #wrap-attach:
        #    ResidueModification
        ARTIFACT
        HYPOTHETICAL
        NATURAL
        POSTTRANSLATIONAL
        MULTIPLE
        CHEMICAL_DERIVATIVE
        ISOTOPIC_LABEL
        PRETRANSLATIONAL
        OTHER_GLYCOSYLATION
        NLINKED_GLYCOSYLATION
        AA_SUBSTITUTION
        OTHER
        NONSTANDARD_RESIDUE
        COTRANSLATIONAL
        OLINKED_GLYCOSYLATION
        UNKNOWN
        NUMBER_OF_SOURCE_CLASSIFICATIONS

