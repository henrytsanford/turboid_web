FROM continuumio/miniconda3

WORKDIR /app

COPY environment.yml .
RUN conda env create -f environment.yml

COPY turboid_notebook.py turboid_differential_expression.csv ./

# Expose port (Cloud Run uses PORT env variable)
EXPOSE 8080

# Run marimo with proper host and port settings
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "marimo", "marimo", "run", "turboid_notebook.py", "--host", "0.0.0.0", "--port", "8080"]