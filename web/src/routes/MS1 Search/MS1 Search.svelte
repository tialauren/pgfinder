<script lang="ts">
  import {
    ProgressBar,
    getModalStore,
    type ModalSettings,
  } from "@skeletonlabs/skeleton";
  import AdvancedOptions from "./AdvancedOptions.svelte";
  import MassDatabaseUploader from "./MassDatabaseUploader.svelte";
  import MsDataUploader from "./MsDataUploader.svelte";
  import ErrorModal from "../ErrorModal.svelte";
  import { CUSTOM_SPECIES_VALUE, defaultPythonState } from "$lib/constants";
  import PGFinder from "$lib/pgfinder.ts?worker";
  import fileDownload from "js-file-download";
  import { onMount } from "svelte";

  export let version: string | undefined;

  // Get the Error Modal Store
  const modalStore = getModalStore();

  // Declare component state
  let state: PythonState = { ...defaultPythonState };
  let allowedModifications: Array<string>;
  let massLibraries: MassLibraryIndex;
  let speciesIndex: SpeciesIndex;
  let loading = true;
  let processing = false;
  let advancedMode = false;
  let previewResult: CustomRulePreviewResult | undefined;
  let previewLoading = false;

  // Start PGFinder and attach callbacks
  let pgfinder: Worker;
  onMount(() => {
    pgfinder = new PGFinder();
    pgfinder.onmessage = ({ data: msg }) => {
      switch (msg.type) {
        case "Ready":
          version = msg.version;
          allowedModifications = msg.allowedModifications;
          massLibraries = msg.massLibraries;
          speciesIndex = msg.speciesIndex;
          loading = false;
          break;
        case "Result":
          fileDownload(msg.blob, msg.filename);
          processing = false;
          break;
        case "Error": {
          const modal: ModalSettings = {
            type: "component",
            component: {
              ref: ErrorModal,
              props: {
                message: msg.message,
              },
            },
          };
          modalStore.trigger(modal);
          processing = false;
          previewLoading = false;
          break;
        }
        case "Preview":
          previewResult = msg;
          previewLoading = false;
          break;
      }
    };
  });

  // Send data to PGFinder for processing
  function runAnalysis() {
    pgfinder.postMessage(state);
    processing = true;
  }

  // Ask PGFinder to show example structures matched by the custom donor/acceptor patterns
  function requestPreview() {
    const req: PGFPreviewReq = { kind: "PreviewCustomRule", state };
    pgfinder.postMessage(req);
    previewLoading = true;
  }

  // Ask PGFinder to generate (and download) just the theoretical dimer list,
  // without matching it against the raw MS1 data
  function generateTheoreticalDimers() {
    const req: PGFGenerateTheoreticalDimersReq = {
      kind: "GenerateTheoreticalDimers",
      state,
    };
    pgfinder.postMessage(req);
    processing = true;
  }

  // A previously-confirmed preview is invalidated as soon as the patterns or
  // bridge constraints change, so a stale "this is fine" can't carry over.
  // eslint-disable-next-line @typescript-eslint/no-unused-vars -- args exist only so `$:` below tracks them as dependencies
  function invalidateCustomRulePreview(...dependencies: unknown[]) {
    previewResult = undefined;
  }
  $: invalidateCustomRulePreview(
    state.customDonorPattern,
    state.customAcceptorPattern,
    state.customAcceptorBridgeType,
    state.customMinGlycineBridge,
    state.customMaxGlycineBridge,
  );

  $: customRuleReady =
    !state.enableDimers ||
    state.species !== CUSTOM_SPECIES_VALUE ||
    (previewResult !== undefined && !previewResult.error);

  // Reactively compute if PGFinder is ready
  $: ready =
    !loading &&
    !processing &&
    state.msData !== undefined &&
    state.massLibrary !== undefined &&
    (!state.enableDimers || state.species !== undefined) &&
    customRuleReady;

  // Theoretical dimers can only be generated once dimer matching is actually configured
  $: theoreticalDimersReady = ready && state.enableDimers;

  // Reactively adapt the UI when entering advanced mode
  $: uiWidth = advancedMode ? "md:w-[40rem]" : "";

  // It's nice to animate the width when opening and closing advanced mode, but
  // it seems like animating the opening leads to some jittery animations, so
  // this is just enables the animation on close. If browsers ever put
  // transitions in their own threads, then maybe this will look nice...
  $: animateWidth = !advancedMode ? "transition-all" : "";
</script>

<div class="flex flex-col items-center">
  <h3 class="pb-1 h3">MS Analysis</h3>
  <div
    class="card m-2 w-[20rem] {uiWidth} max-w-[90%] {animateWidth}"
    data-testid="MS1 Search"
  >
    <section class="flex flex-col space-y-4 p-4">
      <MsDataUploader bind:value={state.msData} />

      <MassDatabaseUploader bind:value={state.massLibrary} {massLibraries} />

      <AdvancedOptions
        bind:enabledModifications={state.enabledModifications}
        bind:ppmTolerance={state.ppmTolerance}
        bind:cleanupWindow={state.cleanupWindow}
        bind:consolidationPpm={state.consolidationPpm}
        bind:enableDimers={state.enableDimers}
        bind:species={state.species}
        bind:permissiveMode={state.permissiveMode}
        bind:customDonorPattern={state.customDonorPattern}
        bind:customAcceptorPattern={state.customAcceptorPattern}
        bind:customLosesTerminalAla={state.customLosesTerminalAla}
        bind:customAcceptorBridgeType={state.customAcceptorBridgeType}
        bind:customMinGlycineBridge={state.customMinGlycineBridge}
        bind:customMaxGlycineBridge={state.customMaxGlycineBridge}
        bind:donorAbundanceThreshold={state.donorAbundanceThreshold}
        bind:advancedMode
        {allowedModifications}
        {speciesIndex}
        {previewResult}
        {previewLoading}
        {requestPreview}
      />

      <button
        type="button"
        class="btn variant-filled"
        on:click={runAnalysis}
        disabled={!ready}
      >
        Run Analysis
      </button>

      {#if state.enableDimers}
        <button
          type="button"
          class="btn variant-ghost"
          on:click={generateTheoreticalDimers}
          disabled={!theoreticalDimersReady}
        >
          Generate Theoretical Dimers
        </button>
      {/if}

      {#if processing}
        <ProgressBar />
      {/if}
    </section>
  </div>
</div>
