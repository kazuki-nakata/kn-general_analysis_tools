from qgis.core import (
    QgsProject,
    QgsColorRampShader,
    QgsSingleBandPseudoColorRenderer,
    QgsRasterShader,
    QgsStyle
)
from PyQt5.QtGui import QColor

# 設定
min_val = 140
max_val = 290
band = 1  # 単一バンドラスタ

# TurboカラーマップをQGISスタイルから取得
style = QgsStyle().defaultStyle()
turbo_ramp = style.colorRamp('Turbo')  # QGIS 3.16以降に含まれている

# Turboカラーマップを等間隔でサンプリングして ColorRampItemList を生成
def generate_color_ramp_items(ramp, n=100):
    items = []
    for i in range(n + 1):
        value = min_val + (max_val - min_val) * i / n
        color = ramp.color(float(i) / n)
        items.append(QgsColorRampShader.ColorRampItem(value, color))
    return items

color_items = generate_color_ramp_items(turbo_ramp)

# 全ラスタレイヤに適用
for layer in QgsProject.instance().mapLayers().values():
    if layer.type() == layer.RasterLayer:
        shader = QgsRasterShader()
        ramp_shader = QgsColorRampShader()
        ramp_shader.setColorRampItemList(color_items)
        ramp_shader.setColorRampType(QgsColorRampShader.Interpolated)
        shader.setRasterShaderFunction(ramp_shader)

        renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), band, shader)
        layer.setRenderer(renderer)
        layer.triggerRepaint()

print("Turboカラーランプ（155–180）を全ラスターレイヤに一括適用しました。")
