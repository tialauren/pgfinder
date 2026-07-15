<script lang="ts">
  import { Accordion, AccordionItem } from "@skeletonlabs/skeleton";
  import ModificationSelector from "./ModificationSelector.svelte";
  import DimerOptions from "./DimerOptions.svelte";
  import Tooltip from "../Tooltip.svelte";
  export let enabledModifications: Array<string>;
  export let allowedModifications: Array<string> | undefined;
  export let ppmTolerance: number;
  export let cleanupWindow: number;
  export let consolidationPpm: number;
  export let advancedMode: boolean;
  export let enableDimers: boolean;
  export let species: string | undefined;
  export let permissiveMode: boolean;
  export let speciesIndex: SpeciesIndex | undefined;
  export let customDonorPattern: string;
  export let customAcceptorPattern: string;
  export let customLosesTerminalAla: boolean;
  export let customAcceptorBridgeType: "none" | "glycine" | "dasp";
  export let customMinGlycineBridge: number;
  export let customMaxGlycineBridge: number;
  export let donorAbundanceThreshold: number;
  export let previewResult: CustomRulePreviewResult | undefined;
  export let previewLoading: boolean;
  export let requestPreview: () => void;
</script>

<Accordion class="w-full">
  <AccordionItem bind:open={advancedMode}>
    <svelte:fragment slot="summary">Advanced Options</svelte:fragment>
    <svelte:fragment slot="content">
      <div class="grid md:grid-cols-2 md:gap-8">
        <ModificationSelector
          bind:value={enabledModifications}
          {allowedModifications}
        />
        <div
          class="flex flex-col justify-between aspect-square overflow-y-auto"
        >
          <div class="flex flex-col items-center">
            <h5 class="pb-1 h5">PPM Tolerance</h5>
            <input
              bind:value={ppmTolerance}
              class="input"
              type="number"
              step="1"
              min="0"
            />
          </div>

          <div class="flex flex-col items-center">
            <h5 class="pb-1 h5">
              Cleanup Window
              <Tooltip style="inline ml-1" type="info">
                Set time window for in-source decay and salt adduct cleanup
              </Tooltip>
            </h5>
            <input
              bind:value={cleanupWindow}
              class="input"
              type="number"
              step="0.1"
              min="0"
            />
          </div>

          <div class="flex flex-col items-center">
            <h5 class="pb-1 h5">
              Consolidation PPM
              <Tooltip style="inline ml-1" type="info">
                During consolidation, structures with the lowest absolute ppm
                are selected over those farther from the theoretical mass.
                However, if two or more matches have a theoretical mass less
                than the consolidation ppm apart, then those matches are
                retained, leaving several possible matches.
              </Tooltip>
            </h5>
            <input
              bind:value={consolidationPpm}
              class="input"
              type="number"
              step="1"
              min="0"
              max={ppmTolerance}
            />
          </div>
        </div>
      </div>

      <DimerOptions
        bind:enableDimers
        bind:species
        bind:permissiveMode
        bind:customDonorPattern
        bind:customAcceptorPattern
        bind:customLosesTerminalAla
        bind:customAcceptorBridgeType
        bind:customMinGlycineBridge
        bind:customMaxGlycineBridge
        bind:donorAbundanceThreshold
        {speciesIndex}
        {previewResult}
        {previewLoading}
        {requestPreview}
      />
    </svelte:fragment>
  </AccordionItem>
</Accordion>
