# MaraTon Challenge 1
## Summary
* The contest aimed to find an efficient way to losslessly compress TON's Bag of Cells (BoC) within a strict time limit of 2 seconds per test case.
* I secured **22nd** place in the contest, consistently holding around the **15th** rank for much of the competition despite a tight schedule during my **_college finals_**.

---

## **My Final Approach**

1. **Compressor**:

   - My final compressor is primarily based on **PAQ8**, merging elements from different versions.
   - The base code was from **tangelo v1.0**, enhanced with models from stronger compressors like **paq8pxd**.

2. **Customizations**:

   - Modified almost all models to better match the BoC data. Some parts of the data, when included in models, made the compression worse, so they were excluded.
   - Fine-tuned constants and parameters through trial and error to optimize performance.

3. **APM Optimization**:

   - Refined contexts and weights in the **Adaptive Probability Maps (APMs)** to better suit the data, achieving improved compression ratios.

4. **Time Limit Fit**:

   - Adjusted compression stress based on input size to fit the 2-second time limit.

---

## **The Journey**

### **Initial Experiments**

- **Huffman Coding**:

  - Preprocessed input files and observed high-frequency bytes even after **lz4** compression.
  - Implemented **Huffman Coding**, which slightly improved compression.

- **Iterative LZ4**:

  - Attempted multiple iterations of lz4 compression but found diminishing returns as the overhead negated any benefits.

- **Algorithm Comparison**:

  - Ran a Python script to test various compression algorithms on the given 25 samples.
  - Results highlighted **GZIP** and **LZMA** as strong contenders, with **LZMA** slightly outperforming.

### **Mid-Contest Improvements**

- **Minimizing Library Size**:

  - Initially used **[miniz-cpp](https://github.com/tfussell/miniz-cpp)**, but its size exceeded the 65KB submission limit.
  - Refined and minified the library, reducing it from 250KB to 65KB—an almost impossible feat!
  - Eventually submitted a solution that improved my rank from **64th to 52nd**.

- **Switch to [TinyLZMA](https://github.com/WangXuan95/TinyLZMA)**:

  - Found **[TinyLZMA](https://github.com/WangXuan95/TinyLZMA)**, which significantly boosted my rank to **23rd** after fine-tuning.
  - With further fine-tuning, I managed to reach the **20th rank**!

- **Introduction to PAQ**:

  - Discovered **Matt Mahoney’s Data Compression website** and explored various PAQ versions.
  - Settled on **tangelo v1.0** for its balance of time and compression ratio. After removing less useful context models and dynamically activating models based on input size, I reached **14th place**.

### **Fine-Tuning and Refinements**

- **Enhancing Models**:

  - Integrated more complex context models from **paq8pxd**, modifying them to fit the BoC data.
  - Adjusted parameters and constants for further improvements.

- **Missed Opportunities**:

  - Realized deeper analysis of the BoC data structure could unlock additional improvements.
  - Unfortunately, my **college finals** limited the time I could dedicate to this.

### **Final Rank**

- Maintained a position in the **top 20** until the last day, finishing in **22nd place** after the final tests.
- The downfall was very fast, and I even suspected cheating at one point, and wrote this [blog](https://codeforces.com/blog/entry/138397), and I want to thank the judges for taking my concerns into consideration, investigating the matter, and clarifying it..

---

## **Failed Trials**

- **Saving Neural Network and Context States**:

  - While beneficial for compression, the resulting model size exceeded the upload limit.

- **Warm-Up Compression**:

  - Pre-trained the neural network using a small sample input, but the improvement was minor.
  - I now feel I didn't give this approach a serious try, as many other leaders mentioned using this technique. A stronger attempt might have yielded better results.

- **BWT Preprocessing**:

  - Preprocessing with **Burrows-Wheeler Transform** worsened the compression ratio.

- **Strong vs. Moderate Compression**:

  - For large files, I tried heavily compressing the initial part of files and lightly compressing the rest but achieved worse results overall.
  - For example, instead of compressing a 300KB file to 200KB, I tried compressing the first 150KB to 75KB and the rest 150KB to 110KB, yielding 185KB. However, this approach was less effective than compressing the whole file moderately.

- **Static Bits in Metadata**:

  - I removed some static bits in the metadata of each cell that were always `0`. However, reading/writing bits instead of bytes introduced significant overhead for very minimal enhancement.

- **Testing Various Algorithms**:

  - Tried almost all promising compression methods mentioned on **[Matt Mahoney’s website](https://mattmahoney.net/dc/)**. While some had weaker compression, others were too slow. **PAQ** was the sweet spot.

- **LLM Compression**:

  - Experimented with **LLM compression** using the **[GMIX library](https://github.com/byronknoll/gmix)**, which offered better performance but was extremely slow. Additionally, the created model could not be uploaded to Codeforces.

---

## **Unimplemented Ideas (No Time)**

- **Splitting Metadata and Payload**:

  - Separate compression for metadata and payload could improve ratios. (Now, after reading some of the other winners' approaches, this idea did have a good impact).

- **Analyzing BoC Data Patterns**:

  - Deeper preprocessing and reorganization of serialized data may yield better compression.

---

## **Advice to me for Future Compression Challenges (Learned the Hard Way)**

1. **Set Ambitious Goals**:

   - If aiming for the **top 20**, work toward a **top 10** solution to account for last-minute submissions.

2. **Understand the Data**:

   - Invest whatever time is needed to understand the underlying data patterns—how to work with it, traverse it, manipulate it, disassemble it, and reassemble it. Become fluent with the data to avoid treating it as a black box.

3. **Prioritize Preprocessing**:

   - Preprocessing data effectively can be more impactful than just optimizing compression.

4. **Serialization Isn't Only Compact**:

   - A serialization being compact doesn’t mean it’s ready for compression. Resorting similar parts together can improve compression. Improving serialization **doesn’t necessarily mean reducing the number of bytes used**!

5. **Examine All Resources**:

   - Carefully examine the resources provided by problem authors, even if they seem irrelevant. Near the end of the competition, I discovered that some code examples were provided in the Mac and Linux versions of the problem (that I didn't download, because I use Windows), these would have been very helpful if I found them earlier.

6. **Search For Easter Eggs xD**:

   - Investigate easter eggs or hidden hints in test cases or the problem statement. Someone got **\$20** from a test case!

---

## **Contest Links**

- [Contest Page](https://codeforces.com/contest/2054)
- [Discussion Blog](https://codeforces.com/blog/entry/137533)

---

Much thanks to [Matt Mahoney’s Data Compression website](https://mattmahoney.net/dc/), which was invaluable throughout this journey.   
Also, a big thank you to the TON team for organizing such an engaging and well-structured contest. It was a fantastic opportunity to push my limits, learn, and grow!
