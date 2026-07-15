<script lang="ts">
  import {
    ListBox,
    ListBoxItem,
    SlideToggle,
    ProgressRadial,
    RadioGroup,
    RadioItem,
  } from "@skeletonlabs/skeleton";
  import Tooltip from "../Tooltip.svelte";
  import { CUSTOM_SPECIES_VALUE } from "$lib/constants";

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

  $: isCustom = species === CUSTOM_SPECIES_VALUE;
</script>

<div class="flex flex-col items-center space-y-2 w-full pt-4">
  <SlideToggle
    name="enable-dimers"
    bind:checked={enableDimers}
    active="bg-primary-500"
  >
    Enable Dimer Matching
  </SlideToggle>

  {#if enableDimers}
    <h5 class="pb-1 h5">Species</h5>
    {#if speciesIndex !== undefined}
      <ListBox class="w-full">
        {#each Object.entries(speciesIndex) as [code, name]}
          <ListBoxItem bind:group={species} name="species" value={code}>
            <i>{name}</i>
          </ListBoxItem>
        {/each}
        <ListBoxItem
          bind:group={species}
          name="species"
          value={CUSTOM_SPECIES_VALUE}
        >
          Custom...
        </ListBoxItem>
      </ListBox>
    {/if}

    <div class="w-full pt-2">
      <h5 class="pb-1 h5">
        Donor Abundance Threshold (%)
        <Tooltip style="inline ml-1" type="info">
          Donors are selected from most to least abundant, structurally eligible
          monomer, until their combined intensity covers this percentage of the
          eligible pool's total intensity. Lower values select fewer, more
          dominant donors; 100% selects every eligible donor.
        </Tooltip>
      </h5>
      <input
        bind:value={donorAbundanceThreshold}
        class="input w-24"
        type="number"
        min="0"
        max="100"
        step="1"
      />
    </div>

    {#if isCustom}
      <div class="flex flex-col w-full space-y-2 pt-2">
        <div>
          <h5 class="pb-1 h5">
            Donor Pattern
            <Tooltip style="inline ml-1" type="info">
              Regex matched against the bridge-stripped stem sequence (e.g.
              "^..JAA$" for a pentapeptide ending D-Ala-D-Ala with mDAP at
              position 3).
            </Tooltip>
          </h5>
          <input
            bind:value={customDonorPattern}
            class="input"
            type="text"
            placeholder="^..JAA$"
          />
        </div>

        <div>
          <h5 class="pb-1 h5">
            Acceptor Pattern
            <Tooltip style="inline ml-1" type="info">
              Regex matched against the bridge-stripped stem sequence.
            </Tooltip>
          </h5>
          <input
            bind:value={customAcceptorPattern}
            class="input"
            type="text"
            placeholder="^..J"
          />
        </div>

        <div>
          <h5 class="pb-1 h5">
            Acceptor Bridge
            <Tooltip style="inline ml-1" type="info">
              Some species crosslink through a peptide bridge attached to the
              acceptor's side chain (e.g. S. aureus's pentaglycine bridge, or E.
              faecalis's single D-Asp/D-Asn bridge) rather than bonding
              directly. Choose "None" for a direct, bridge-less crosslink.
            </Tooltip>
          </h5>
          <RadioGroup>
            <RadioItem
              bind:group={customAcceptorBridgeType}
              name="bridge-type"
              value="none"
            >
              None
            </RadioItem>
            <RadioItem
              bind:group={customAcceptorBridgeType}
              name="bridge-type"
              value="glycine"
            >
              Glycine
            </RadioItem>
            <RadioItem
              bind:group={customAcceptorBridgeType}
              name="bridge-type"
              value="dasp"
            >
              D-Asp
            </RadioItem>
          </RadioGroup>

          {#if customAcceptorBridgeType === "glycine"}
            <div class="flex items-center justify-center space-x-2 pt-2">
              <span>Bridge length</span>
              <input
                bind:value={customMinGlycineBridge}
                class="input w-16"
                type="number"
                min="1"
                max={customMaxGlycineBridge}
              />
              <span>to</span>
              <input
                bind:value={customMaxGlycineBridge}
                class="input w-16"
                type="number"
                min={customMinGlycineBridge}
              />
              <span>glycines</span>
            </div>
          {/if}
        </div>

        <label class="flex items-center space-x-2">
          <SlideToggle
            name="custom-loses-terminal-ala"
            bind:checked={customLosesTerminalAla}
            size="sm"
            active="bg-primary-500"
          />
          <span>
            Donor loses terminal D-Ala
            <Tooltip style="inline ml-1" type="info">
              Whether the donor loses a terminal D-Ala residue during crosslink
              formation, as in the standard 4-3 mechanism.
            </Tooltip>
          </span>
        </label>

        <button
          type="button"
          class="btn variant-filled"
          on:click={requestPreview}
          disabled={!customDonorPattern ||
            !customAcceptorPattern ||
            previewLoading}
        >
          {#if previewLoading}
            <ProgressRadial width="w-4" />
          {:else}
            Preview Matches
          {/if}
        </button>

        {#if previewResult !== undefined}
          {#if previewResult.error}
            <p class="text-error-500">{previewResult.error}</p>
          {:else}
            <div class="text-sm">
              <p>
                <b>Donor pattern</b> matches {previewResult.donorCount} structure{previewResult.donorCount ===
                1
                  ? ""
                  : "s"}
              </p>
              {#each previewResult.donorMatches as match}
                <p class="pl-2 truncate">{match}</p>
              {/each}

              <p class="pt-2">
                <b>Acceptor pattern</b> matches {previewResult.acceptorCount} structure{previewResult.acceptorCount ===
                1
                  ? ""
                  : "s"}
              </p>
              {#each previewResult.acceptorMatches as match}
                <p class="pl-2 truncate">{match}</p>
              {/each}
            </div>
          {/if}
        {/if}
      </div>
    {:else}
      <div class="flex items-center space-x-2 pt-2">
        <SlideToggle
          name="permissive-mode"
          bind:checked={permissiveMode}
          size="sm"
          active="bg-primary-500"
        />
        <span>
          Permissive Mode
          <Tooltip style="inline ml-1" type="info">
            Allows bridge-length variants for novel crosslink discovery, instead
            of enforcing strict literature rules.
          </Tooltip>
        </span>
      </div>
    {/if}
  {/if}
</div>
