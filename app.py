# app.py
import streamlit as st
import pandas as pd
import altair as alt
from alphagenome.data import gene_annotation, genome, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

# --- Page Configuration ---
st.set_page_config(
    page_title="AlphaGenome Variant Visualizer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- Global Variables & Constants ---
ORGANISM_MAP = {
    'human': dna_client.Organism.HOMO_SAPIENS,
    'mouse': dna_client.Organism.MUS_MUSCULUS,
}

HG38_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)
MM10_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'mm10/gencode.vM23.annotation.gtf.gz.feather'
)

# --- Caching Functions ---

@st.cache_resource
def get_dna_model():
    """
    Creates and caches the DNA model client using the API key from Streamlit secrets.
    """
    api_key = st.secrets.get("ALPHAGENOME_API_KEY")
    if not api_key:
        st.error("AlphaGenome API key not found. Please add it to your Streamlit secrets.")
        st.info("For local testing, create a file at .streamlit/secrets.toml with the content: \n\nALPHAGENOME_API_KEY = \"YOUR_API_KEY_HERE\"")
        st.stop()
    try:
        return dna_client.create(api_key)
    except Exception as e:
        st.error(f"Failed to create DNA model client. Please check your API key. Error: {e}")
        st.stop()

@st.cache_data
def get_transcript_extractors(organism_str):
    """
    Loads and caches gene annotation data and creates transcript extractors.
    """
    organism = ORGANISM_MAP[organism_str]
    with st.spinner('Loading gene annotation...'):
        if organism == dna_client.Organism.HOMO_SAPIENS:
            gtf_path = HG38_GTF_FEATHER
        elif organism == dna_client.Organism.MUS_MUSCULUS:
            gtf_path = MM10_GTF_FEATHER
        else:
            st.error(f'Unsupported organism: {organism}')
            st.stop()

        gtf = pd.read_feather(gtf_path)
        gtf_transcript = gene_annotation.filter_transcript_support_level(
            gene_annotation.filter_protein_coding(gtf), ['1']
        )
        transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)
        gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(
            gtf_transcript
        )
        longest_transcript_extractor = transcript.TranscriptExtractor(
            gtf_longest_transcript
        )
        return transcript_extractor, longest_transcript_extractor

@st.cache_data
def predict_variant_cached(_dna_model, interval_str, variant_str, organism_str):
    """
    Cache wrapper for dna_model.predict_variant to avoid re-running predictions.
    The _dna_model argument is used to invalidate the cache if the model changes,
    but it's not directly used in the function body.
    """
    interval = genome.Interval.from_str(interval_str)
    variant = genome.Variant.from_str(variant_str)
    organism = dna_client.Organism[organism_str]
    
    # This is the correct prediction function based on the working Colab notebook
    return _dna_model.predict_variant(
      interval=interval,
      variant=variant,
      organism=organism,
      requested_outputs=[*dna_client.OutputType],
      ontology_terms=['EFO:0001187', 'EFO:0002067', 'EFO:0002784'],
  )

# --- Sidebar UI ---
st.sidebar.title('Variant Visualizer')
st.sidebar.header('Variant and Plotting Options')

organism_choice = st.sidebar.selectbox('Organism', ['human', 'mouse'])
variant_chromosome = st.sidebar.text_input('Chromosome', 'chr22')
variant_position = st.sidebar.number_input('Position', value=36201698, format="%d")
variant_reference_bases = st.sidebar.text_input('Reference Bases', 'A')
variant_alternate_bases = st.sidebar.text_input('Alternate Bases', 'C')
sequence_length_choice = st.sidebar.selectbox('Sequence Length', ["2KB", "16KB", "100KB", "500KB", "1MB"], index=4)

st.sidebar.header('Plotting Options')
st.sidebar.subheader('Output Types')
plot_rna_seq = st.sidebar.checkbox('RNA-SEQ', value=True)
plot_splice_sites = st.sidebar.checkbox('SPLICE_SITES', value=True)
plot_splice_junctions = st.sidebar.checkbox('SPLICE_JUNCTIONS', value=True)
plot_cage = st.sidebar.checkbox('CAGE', value=True)
plot_atac = st.sidebar.checkbox('ATAC', value=False)
plot_dnase = st.sidebar.checkbox('DNASE', value=False)
plot_chip_histone = st.sidebar.checkbox('CHIP_HISTONE', value=False)
plot_chip_tf = st.sidebar.checkbox('CHIP_TF', value=False)
plot_contact_maps = st.sidebar.checkbox('CONTACT_MAPS', value=False)

st.sidebar.subheader('Gene Annotation')
plot_gene_annotation = st.sidebar.checkbox('Plot Gene Annotation', value=True)
plot_longest_transcript_only = st.sidebar.checkbox('Plot Longest Transcript Only', value=True)

st.sidebar.subheader('DNA Strand Filter')
strand_filter = st.sidebar.radio('Filter to Strand', ('None', 'Positive', 'Negative'), index=0)

st.sidebar.subheader('Visualization Settings')
plot_interval_width = st.sidebar.slider('Plot Interval Width (bp)', min_value=1024, max_value=200000, step=1024, value=43008)
plot_interval_shift = st.sidebar.slider('Plot Interval Shift (bp)', min_value=-524288, max_value=524288, step=2048, value=0)
ref_color = st.sidebar.color_picker('Reference Color', value='#808080')
alt_color = st.sidebar.color_picker('Alternate Color', value='#FF0000')
ref_alt_colors = {'REF': ref_color, 'ALT': alt_color}

# --- Main App Logic ---
dna_model = get_dna_model()

if st.sidebar.button('Visualize Variant', use_container_width=True):
    st.header(f"Visualizing variant: {variant_chromosome}:{variant_position}:{variant_reference_bases}>{variant_alternate_bases}")

    variant = genome.Variant(
        chromosome=variant_chromosome,
        position=int(variant_position),
        reference_bases=variant_reference_bases,
        alternate_bases=variant_alternate_bases,
    )

    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f'SEQUENCE_LENGTH_{sequence_length_choice}']
    interval = variant.reference_interval.resize(sequence_length)

    transcript_extractor, longest_transcript_extractor = get_transcript_extractors(organism_choice)

    output = predict_variant_cached(
        dna_model,
        interval_str=str(interval),
        variant_str=str(variant),
        organism_str=ORGANISM_MAP[organism_choice].name,
    )

    ref, alt = output.reference, output.alternate

    if strand_filter == 'Positive':
        ref = ref.filter_to_strand(strand='+')
        alt = alt.filter_to_strand(strand='+')
    elif strand_filter == 'Negative':
        ref = ref.filter_to_strand(strand='-')
        alt = alt.filter_to_strand(strand='-')

    components = []
    if plot_gene_annotation:
        extractor = longest_transcript_extractor if plot_longest_transcript_only else transcript_extractor
        transcripts = extractor.extract(interval)
        components.append(plot_components.TranscriptAnnotation(transcripts))

    plot_map = {
        plot_atac: (ref.atac, alt.atac, 'ATAC'),
        plot_cage: (ref.cage, alt.cage, 'CAGE'),
        plot_chip_histone: (ref.chip_histone, alt.chip_histone, 'CHIP_HISTONE'),
        plot_chip_tf: (ref.chip_tf, alt.chip_tf, 'CHIP_TF'),
        plot_contact_maps: (ref.contact_maps, alt.contact_maps, 'CONTACT_MAPS'),
        plot_dnase: (ref.dnase, alt.dnase, 'DNASE'),
        plot_rna_seq: (ref.rna_seq, alt.rna_seq, 'RNA_SEQ'),
        plot_splice_junctions: (ref.splice_junctions, alt.splice_junctions, 'SPLICE_JUNCTIONS'),
        plot_splice_sites: (ref.splice_sites, alt.splice_sites, 'SPLICE_SITES'),
    }

    for should_plot, (ref_data, alt_data, output_type) in plot_map.items():
        if should_plot:
            if ref_data is None or ref_data.values.shape[-1] == 0:
                st.warning(f'No tracks exist for {output_type} with the current filters.')
                continue

            if output_type == 'CONTACT_MAPS':
                components.append(plot_components.ContactMapsDiff(tdata=alt_data - ref_data))
            elif output_type == 'SPLICE_JUNCTIONS':
                components.append(plot_components.Sashimi(ref_data, ylabel_template='REF: {track}'))
                components.append(plot_components.Sashimi(alt_data, ylabel_template='ALT: {track}'))
            else:
                components.append(plot_components.OverlaidTracks(
                    tdata={'REF': ref_data, 'ALT': alt_data},
                    colors=ref_alt_colors,
                    ylabel_template='{track}'
                ))

    if not components:
        st.warning("No data to plot. Please select at least one output type to visualize.")
    elif plot_interval_width > interval.width:
        st.error(f'Plot Interval Width ({plot_interval_width}) must be less than Sequence Length ({interval.width}).')
    else:
        with st.spinner('Generating plot...'):
            plot = plot_components.plot(
                components=components,
                interval=interval.shift(plot_interval_shift).resize(plot_interval_width),
                annotations=[plot_components.VariantAnnotation([variant])],
            )
            st.altair_chart(plot, use_container_width=True)
else:
    st.info("Configure your variant in the sidebar and click 'Visualize Variant' to begin.")

