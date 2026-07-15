export const defaultVersions = {
  WebUI: WEBUI_VERSION,
  Smithereens: undefined,
  PGFinder: undefined,
};

export const defaultPythonState = {
  msData: undefined,
  massLibrary: undefined,
  enabledModifications: [],
  ppmTolerance: 10,
  cleanupWindow: 0.5,
  consolidationPpm: 1,
  enableDimers: false,
  species: undefined,
  permissiveMode: false,
  customDonorPattern: "",
  customAcceptorPattern: "",
  customLosesTerminalAla: false,
  customAcceptorBridgeType: "none" as const,
  customMinGlycineBridge: 5,
  customMaxGlycineBridge: 5,
  donorAbundanceThreshold: 90,
};

// Must match CUSTOM_SPECIES_VALUE in lib/pgfinder/gui/shim.py
export const CUSTOM_SPECIES_VALUE = "__custom__";
