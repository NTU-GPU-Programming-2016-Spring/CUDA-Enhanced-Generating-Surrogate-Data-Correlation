## 檔案路徑說明

林教授提供的五項資源，解壓後如下：

+ fsaverage/ (此處不提供)
+ lol_stc/ (此處不提供)
+ matlab_functions/
+ render_example/ (此處不提供)
+ subject_folder.xlsx

matlab_functions 會新增程式，並於下方補充說明。基本上，這裡提供的檔案，全部可直接覆蓋教授所提供的檔案。

## 環境建置

以我來說，上述的檔案的根目錄為：

`C:\Salmon_Projects\GPU_final`

於 Matlab 輸入，增加 Matlab **Path**：

`path('C:\Salmon_Projects\GPU_final\matlab_functions', path)`

於 Matlab 輸入，檢查 **Path**：

`path()`

於 Matlab 輸入，增加環境變數 **SUBJECTS_DIR**：

`setenv('SUBJECTS_DIR', 'C:\Salmon_Projects\GPU_final')`

於 Matlab 輸入，檢查環境變數 **SUBJECTS_DIR**：

`getenv('SUBJECTS_DIR')`

## 相關操作

### 以 Matlab 讀取腦半球活動紀錄

使用 `inverse_read_stc.m` 的函數，以 `s01_120415_2_fsaverage_fmcstcs-lh.stc` 為例，return 為一 martix 物件。

```
inverse_read_stc('../lol_stc/s01_120415/unpack/bold/002/s01_120415_2_fsaverage_fmcstcs-lh.stc')
```

---

### .stc 格式轉 .csv

檔案名稱: `convert2csv.m`。

因為 Matlab column based 的關係，原始資料之 row x column = **時間(秒數) x 大腦活動數據** = **約 450 x 10242**。以下提供了兩種轉換方式。執行完畢之後，`csv` 會存放在 `lol_stc_csv` 之目錄下。

#### 使用方式

```
% column based
convert2csv('../lol_stc', false);

% row based
convert2csv('../lol_stc', true);
```

#### 輸出目錄結構說明

> ./lol\_stc\_csv/<受試者編號(前)>-<受試時間>-<左右半腦>.csv

e.g.

```
./lol_stc_csv/s01-002-lh.csv
./lol_stc_csv/s01-002-rh.csv
...
```

(written by Salmon 05.11 06')

---

### 輸出全部 stc 大腦圖片

檔案名稱: `render_all_brains.m`。

這一部分會自動 mapping 受試時觀看腳色視角(AD1, SUP1, AD2 and SUP2)。

1st 參數輸入 stc 的 root 資料夾位置；2nd 參數輸入 subject_folder.xlsx 的路徑。執行完畢之後，輸出的圖片 `.png` 會存放在 `lol_stc_brain_image_render` 之目錄下。

#### 使用方式

```
render_all_brains('../lol_stc', '../subject_folder.xlsx');
```

#### 輸出目錄結構說明

> ./lol\_stc_brain\_image\_render/<受試者編號>/<受試時間>-<視角>-<左右半腦>/n.png

e.g.

```
./lol_stc_brain_image_render/s01_120415/002_AD1_lh/1.png
./lol_stc_brain_image_render/s01_120415/002_AD1_lh/2.png
./lol_stc_brain_image_render/s01_120415/006_AD2_rh/1.png
...
```

p.s. 1 暫時未作截頭截尾的部分。

(written by Salmon 05.19 06')

---
