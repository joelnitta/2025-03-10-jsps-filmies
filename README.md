# slides_template

A template repo for rendering presentation slides with [Quarto](https://quarto.org/)

## Rendering slides

First, initialize the repository for publishing:

```
quarto publish gh-pages
```

This creates the `gh-pages` branch and pushes it to GitHub.

Any further commits to `main` will trigger an automatic rendering of the website when pushed to the remote via the [publish.yml](.github/workflows/publish.yml) GitHub Action.

The rendered slides are at <https://joelnitta.github.io/slides_template>

## License

slides_template (c) by Joel H. Nitta
 
slides_template is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

