# app.py
import streamlit as st
import pandas as pd
import altair as alt
import numpy as np
from alphagenome.data import gene_annotation, genome, transcript
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components
import pkg_resources # Import the library to check package versions

# --- Page Configuration ---
st.set_page_config(
    page_title="AlphaGenome Variant Scorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- Debugging Information ---
try:
    ag_version = pkg_resources.get_distribution("alphagenome").version
    st.sidebar.info(f"AlphaGenome Version: `{ag_version}`")
except pkg_resources.DistributionNotFound:
    st.sidebar.error("AlphaGenome library not found.")


# --- Caching Functions ---

@st.cache_resource
def get_dna_model():
    """
    Creates and caches the DNA model client using the API key from Streamlit secrets.
    """
    api_key = st.secrets.get("ALPHAGENOME_API_KEY")
    if not api_key:
        st.error("AlphaGenome API key not found. Please add it to your Streamlit secrets.")
        st.info("Create a file .streamlit/secrets.toml with ALPHAGENOME_API_KEY = 'your_key_here' for local testing, or add it in the advanced settings on Streamlit Community Cloud.")
        st.stop()
    try:
        return dna_client.create(api_key)
    except Exception as e:
        st.error(f"Failed to create DNA model client. Please check your API key. Error: {e}")
        st.stop()

@st.cache_data
def get_gene_annotations(organism_enum):
    """
    Caches the gene annotation data for a given organism.
    NOTE: Using `GeneAnnotationDb` which is expected in newer alphagenome versions.
    """
    try:
        if organism_enum == dna_client.Organism.HOMO_SAPIENS:
            gtf_path = (
                'https://storage.googleapis.com/alphagenome/reference/gencode/'
                'hg38/gencode.v46.annotation.gtf.gz.feather'
            )
        else:
            gtf_path = (
                'https://storage.googleapis.com/alphagenome/reference/gencode/'
                'mm10/gencode.vM23.annotation.gtf.gz.feather'
            )
        
        # Use the newer class name expected in alphagenome>=0.1.0
        return gene_annotation.GeneAnnotationDb(gtf_path)

    except Exception as e:
        st.error(f"Failed to load gene annotation data. Error: {e}")
        return None

# --- Main App UI ---

st.title("🧬 AlphaGenome Variant Scoring & Visualization")
st.write(
    "This app uses the AlphaGenome model to predict the effects of genetic variants "
    "on various genomic modalities and visualize the results."
)

# --- Sidebar for User Inputs ---

with st.sidebar:
    st.header("1. Variant Scoring")

    organism_map = {
        'Human': dna_client.Organism.HOMO_SAPIENS,
        'Mouse': dna_client.Organism.MUS_MUSCULUS,
    }
    organism_label = st.selectbox("Organism", options=list(organism_map.keys()))
    organism = organism_map[organism_label]

    st.subheader("Variant Details")
    col1, col2 = st.columns(2)
    with col1:
        variant_chromosome = st.text_input("Chromosome", "chr22")
    with col2:
        variant_position = st.number_input("Position", value=36201698, step=1, format="%d")

    col3, col4 = st.columns(2)
    with col3:
        variant_reference = st.text_input("Reference", "A")
    with col4:
        variant_alternate = st.text_input("Alternate", "C")

    sequence_length_map = {
        '1MB': dna_client.SEQUENCE_LENGTH_1MB,
        '500KB': dna_client.SEQUENCE_LENGTH_500KB,
        '100KB': dna_client.SEQUENCE_LENGTH_100KB,
        '16KB': dna_client.SEQUENCE_LENGTH_16KB,
        '2KB': dna_client.SEQUENCE_LENGTH_2KB,
    }
    sequence_length_label = st.selectbox("Sequence Length", options=list(sequence_length_map.keys()))
    sequence_length = sequence_length_map[sequence_length_label]

    score_button = st.button("Score Variant", use_container_width=True)

# --- Main Content Area ---

# Initialize session state
if 'scores_df' not in st.session_state:
    st.session_state.scores_df = None
if 'variant' not in st.session_state:
    st.session_state.variant = None
if 'interval' not in st.session_state:
    st.session_state.interval = None

# Check for API key before proceeding
if 'ALPHAGENOME_API_KEY' not in st.secrets:
    st.warning("Please configure your AlphaGenome API key in Streamlit secrets to use this app.")
else:
    if score_button:
        dna_model = get_dna_model()
        try:
            variant = genome.Variant(
                chromosome=variant_chromosome,
                position=variant_position,
                reference_bases=variant_reference,
                alternate_bases=variant_alternate,
            )
            interval = variant.reference_interval.resize(sequence_length)

            with st.spinner("Scoring variant... This may take a moment."):
                variant_scores = dna_model.score_variant(
                    interval=interval,
                    variant=variant,
                    variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
                )
                scores_df = variant_scorers.tidy_scores(variant_scores)

            st.session_state.scores_df = scores_df
            st.session_state.variant = variant
            st.session_state.interval = interval
            st.session_state.organism = organism

        except Exception as e:
            st.error(f"An error occurred during scoring: {e}")
            st.session_state.scores_df = None


    if st.session_state.scores_df is not None:
        st.header("📊 Scoring Results")
        st.info(f"Displaying scores for variant: {st.session_state.variant}")

        st.dataframe(st.session_state.scores_df)

        csv = st.session_state.scores_df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="Download Scores as CSV",
            data=csv,
            file_name=f'{st.session_state.variant}_scores.csv',
            mime='text/csv',
        )

        # --- Visualization Section ---
        with st.sidebar:
            st.header("2. Variant Visualization")
            st.write("Select a track from the results to visualize its predicted effect.")

            unique_tracks = st.session_state.scores_df['track_name'].dropna().unique()
            track_name = st.selectbox("Track to Visualize", options=unique_tracks)

            st.subheader("Plotting Options")
            plot_interval_width = st.number_input("Plot Interval Width (bp)", min_value=100, max_value=10001, value=2001, step=100)
            plot_interval_shift = st.number_input("Plot Interval Shift (bp)", value=0, step=50)

            col5, col6 = st.columns(2)
            with col5:
                separate_strands = st.checkbox("Separate Strands", value=True)
            with col6:
                sashimi_style = st.checkbox("Sashimi Style", value=True)

            visualize_button = st.button("Visualize Variant Effect", use_container_width=True)

        if visualize_button and track_name:
            try:
                with st.spinner("Generating visualization..."):
                    dna_model = get_dna_model()
                    gene_annotations = get_gene_annotations(st.session_state.organism)

                    # Use the newer 'predict_on_batch' method expected in alphagenome>=0.1.0
                    predictions = dna_model.predict_on_batch(
                        inputs={
                            'interval': np.array([str(st.session_state.interval)]),
                            'variant': np.array([str(st.session_state.variant)]),
                        },
                        tracks=[track_name]
                    )
                    # Extract the first (and only) result from the batch
                    ref_data = pd.DataFrame(predictions['ref'][track_name][0])
                    alt_data = pd.DataFrame(predictions['alt'][track_name][0])


                    # Setup plot components
                    components = []
                    
                    if gene_annotations:
                        transcript_extractor = transcript.TranscriptExtractor(gene_annotations)
                        transcripts = transcript_extractor.extract_transcripts(st.session_state.interval)
                        components.append(plot_components.GeneAnnotation(transcripts))
                    
                    ref_alt_colors = {'REF': '#1f77b4', 'ALT': '#ff7f0e'}
                    ylabel_template = '{track}'
                    is_splicing_track = 'SPLICE_SITE_USAGE' in track_name
                    
                    if sashimi_style and is_splicing_track:
                        ref_plot = plot_components.Sashimi(ref_data, ylabel_template='REF: ' + ylabel_template)
                        alt_plot = plot_components.Sashimi(alt_data, ylabel_template='ALT: ' + ylabel_template)
                        components.extend([ref_plot, alt_plot])
                    else:
                        component = plot_components.OverlaidTracks(
                            tdata={'REF': ref_data, 'ALT': alt_data},
                            colors=ref_alt_colors,
                            ylabel_template=ylabel_template,
                            separate_strands=separate_strands,
                        )
                        components.append(component)

                    plot = plot_components.plot(
                        components=components,
                        interval=st.session_state.interval.shift(plot_interval_shift).resize(plot_interval_width),
                        annotations=[plot_components.VariantAnnotation([st.session_state.variant])],
                    )
                    
                    st.header(f"🔬 Visualization for {track_name}")
                    st.altair_chart(plot, use_container_width=True)

            except Exception as e:
                st.error(f"An error occurred during visualization: {e}")
