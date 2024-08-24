import subprocess
import os
import argparse

def generate_embeddings(input_file):
    input_dir = os.path.dirname(os.path.abspath(input_file))
    print("------>",input_dir)
    input_filename = os.path.basename(input_file)
    
    output_file = os.path.join(input_dir, "EmbeddingsData.csv")

    cmd = [
        "docker", "run", "--rm",
        "-v", f"{input_dir}:/input", 
        "scgpt_human",  
        "--input", f"/input/{input_filename}", 
        "--model_dir", "/app" 
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"Embeddings generated successfully. Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error generating embeddings: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Generate gene embeddings using scGPT model')
    parser.add_argument('--input', required=True, help='Path to input expression data CSV file')
    
    args = parser.parse_args()
    
    generate_embeddings(args.input)

if __name__ == "__main__":
    main()
