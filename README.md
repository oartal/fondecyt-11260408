# Sitio web del proyecto Fondecyt de Iniciación 11260408

Sitio estático (HTML + CSS, sin compilación) para GitHub Pages.

## Estructura

```
index.html                     Proyecto: resumen, sociedad, impacto, equipo
estilo.css                     Estilo común
laboratorio/index.html         Presentación de los módulos
laboratorio/remolinos-y-ondas.html   ← reemplazar por el laboratorio real
difusion/index.html            Video, mascotas, actividades, material, notas
difusion/material/             ← poner aquí triptico.pdf, poster.pdf, marcapaginas.pdf
img/                           ← logo-anid.png, logo-ufro.png, imágenes de mascotas
```

## Publicar en GitHub Pages

1. Crear un repositorio público en GitHub (por ejemplo `desertores`).
2. Subir todo el contenido de esta carpeta (botón *Add file → Upload files*).
3. Ir a *Settings → Pages*, elegir *Deploy from a branch*, rama `main`, carpeta `/ (root)`, y guardar.
4. En uno o dos minutos el sitio queda en `https://<usuario>.github.io/desertores/`.

## Pendientes

- Copiar `remolinos-y-ondas.html` desde MaterialDifusion a `laboratorio/`.
- Agregar logos oficiales en `img/` (los del pie se ocultan solos si faltan).
- Video: subirlo a YouTube y pegar el iframe en `difusion/index.html` (GitHub no acepta archivos sobre 100 MB).
- Imágenes de las mascotas y PDFs del material impreso.
