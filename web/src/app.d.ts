declare const WEBUI_VERSION: string;

// FIXME: Remove this! I can just use the native JS `File`!
declare type VirtFile = { name: string; content: ArrayBuffer };

declare type Versions = {
  WebUI: string;
  Smithereens: string | undefined;
  PGFinder: string | undefined;
};

declare type PythonState = {
  msData: Array<VirtFile> | undefined;
  massLibrary: VirtFile | undefined;
  enabledModifications: Array<string>;
  ppmTolerance: number;
  cleanupWindow: number;
  consolidationPpm: number;
  enableDimers: boolean;
  species: string | undefined;
  permissiveMode: boolean;
  customDonorPattern: string;
  customAcceptorPattern: string;
  customLosesTerminalAla: boolean;
  customAcceptorBridgeType: "none" | "glycine" | "dasp";
  customMinGlycineBridge: number;
  customMaxGlycineBridge: number;
  donorAbundanceThreshold: number;
};

declare type CustomRulePreviewResult = {
  donorMatches: Array<string>;
  donorCount: number;
  acceptorMatches: Array<string>;
  acceptorCount: number;
  error: string | undefined;
};

declare type MassLibraryIndex = {
  [index: string]: {
    [index: string]: {
      file: string;
      description: string;
    };
  };
};

declare type SpeciesIndex = {
  [code: string]: string;
};

declare type StructuresIndex = {
  [species: string]: {
    file: string;
    citation: string | undefined;
  };
};

// PGFinder Worker Message Types ===============================================

declare type PGFinderMsg =
  | PGFReadyMsg
  | PGFResultMsg
  | PGFErrorMsg
  | PGFPreviewMsg;

declare type PGFReadyMsg = {
  type: "Ready";
  version: string;
  allowedModifications: Array<string>;
  massLibraries: MassLibraryIndex;
  speciesIndex: SpeciesIndex;
};

declare type PGFResultMsg = {
  type: "Result";
  filename: string;
  blob: Blob;
};

declare type PGFErrorMsg = {
  type: "Error";
  message: string;
};

declare type PGFPreviewMsg = {
  type: "Preview";
} & CustomRulePreviewResult;

declare type PGFPreviewReq = {
  kind: "PreviewCustomRule";
  state: PythonState;
};

declare type PGFGenerateTheoreticalDimersReq = {
  kind: "GenerateTheoreticalDimers";
  state: PythonState;
};

// Smithereens Worker Message Types ============================================

declare type SmithereensReq =
  | SMassReq
  | SMassesReq
  | SValidateReq
  | SFragmentReq
  | SFragmentsReq;

declare type SMassReq = {
  type: "MassReq";
  structure: string;
};

declare type SMassesReq = {
  type: "MassesReq";
  structures: File;
};

declare type SValidateReq = {
  type: "ValidateReq";
  structure: string;
};

declare type SFragmentReq = {
  type: "FragmentReq";
  structure: string;
};

declare type SFragmentsReq = {
  type: "FragmentsReq";
  structures: File;
};

// -----------------------------------------------------------------------------

declare type SmithereensRes =
  | SReady
  | SSingleErr
  | SBulkErr
  | SMassRes
  | SMassesRes
  | SValidateRes
  | SFragmentRes
  | SFragmentsRes;

declare type SReady = {
  type: "Ready";
  version: string;
};

declare type SSingleErr = {
  type: "SingleErr";
};

declare type SBulkErr = {
  type: "BulkErr";
  structure: string;
  line: number;
};

declare type SMassRes = {
  type: "MassRes";
  mass: string;
  smiles: string;
};

declare type SMassesRes = {
  type: "MassesRes";
  filename: string;
  blob: Blob;
};

declare type SValidateRes = {
  type: "ValidateRes";
};

declare type SFragmentRes = {
  type: "FragmentRes";
  filename: string;
  blob: Blob;
};

declare type SFragmentsRes = {
  type: "FragmentsRes";
  filename: string;
  blob: Blob;
};
